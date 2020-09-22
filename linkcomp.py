#---------------------------------------------------------------------
#PACKAGE IMPORTS
#---------------------------------------------------------------------

import sys
import os
import ccdc
from ccdc import search, io, molecule
from ccdc.io import MoleculeReader, CrystalReader, EntryReader, CrystalWriter
from ccdc.descriptors import MolecularDescriptors as MD
from ccdc.descriptors import CrystalDescriptors as CD
from ccdc.crystal import PackingSimilarity
similarity_engine = PackingSimilarity()
import numpy as np
import gemmi
from gemmi import cif
import sympy as sp
import csv

#---------------------------------------------------------------------
#SYSTEM PARAMETERS
#---------------------------------------------------------------------

# Check parameters
if (len(sys.argv) != 3):
    print("Must have 2 parameters: input dir and an output dir!")
    sys.exit()
inputdir = sys.argv[1]
outputdir = sys.argv[2]
# Check they are directories
if not(os.path.isdir(inputdir) & os.path.isdir(outputdir)) :
    print("Parameters are not directories!")
    sys.exit()

#---------------------------------------------------------------------
#CODE FOR CALCULATION OF LINKING NUMBER
#---------------------------------------------------------------------
s,t = sp.symbols('s t')


# triple product function
def trp (x,y,z):
    return np.dot(x,np.cross(y,z))

# denominator functions
def cyc (x,y,z):
    return np.dot(x,y)*np.linalg.norm(z)

#General AT function
def atan(a,b,alpha,d):

    if alpha == 0 or alpha == np.pi or d == 0:
        return 0

    else:
        t = a*b*np.sin(alpha) + d**2*(1/np.tan(alpha))
        m = d*np.sqrt(a**2 + b**2 - (2*a*b*np.cos(alpha))+d**2)

    return np.arctan(t/m)

#Construct invariants of two line segments from endpoints
class Segquad:

    def __init__(self,parameter):

        self.lonestart = np.array(parameter[0][0]).astype(float)
        self.loneend = np.array(parameter[0][1]).astype(float)
        self.ltwostart = np.array(parameter[1][0]).astype(float)
        self.ltwoend = np.array(parameter[1][1]).astype(float)

        #Parameterisation vctors
        self.lonevec = self.loneend - self.lonestart
        self.ltwovec = self.ltwoend - self.ltwostart
        self.startdist = self.ltwostart - self.lonestart

        #Parameterised line equations for segments
        self.loneeq = self.lonestart + t*self.lonevec
        self.ltwoeq = self.ltwostart + s*self.ltwovec

        #Length and cross product for segments
        self.cross = np.cross(self.lonevec, self.ltwovec)
        self.seg1len = np.linalg.norm( self.lonevec )
        self.seg2len = np.linalg.norm( self.ltwovec )

        #symbolic form of gauss integral function
        self.difference = (self.lonestart + t * self.lonevec) - (self.ltwostart + s * self.ltwovec)
        self.distance = sp.sqrt(self.difference[0]**2 + self.difference[1]**2 + self.difference[2]**2)
        self.gaussfunc = np.dot( np.cross( self.lonevec, self.ltwovec ), self.difference) / self.distance**3

        #Angle between segments
        self.angle = np.arccos( np.dot(self.lonevec,self.ltwovec) / (self.seg1len * self.seg2len) )
        #print( 'angle = ' + repr( self.angle ) )
        self.alpha = self.angle/2


    #Calculates the signed distance of a line pair
    def segdist(self):

        return trp(self.lonevec,self.ltwovec,self.startdist)/np.linalg.norm(self.cross)

    #Calculate the a1, a2 co-ordinates of a line pair
    def acoord(self):

        b = self.startdist/((np.sin(self.angle))**2)
        x = self.lonevec/self.seg1len
        y = self.ltwovec/self.seg2len

        a1 = np.dot(((y*np.cos(self.angle))-x),b)
        a2 = np.dot(((y - x*np.cos(self.angle))),b)

        return [a1,a2]

    #Calculates the lk funtion from four separate AT functions
    def segatan(self):

        a1 = self.acoord()[0]
        a2 = self.acoord()[1]
        b1 = self.acoord()[0] + self.seg1len
        b2 = self.acoord()[1] + self.seg2len
        alpha = self.angle
        d = self.segdist()

        return (atan(a1,b2,d,alpha) + atan(a2,b1,d,alpha) - atan(a1,a2,d,alpha) - atan(b1,b2,d,alpha))/(4*np.pi)


#---------------------------------------------------------------------
#GIVEN A LIST OF POINTS p_0, ..., p_n, CALCULATE THE LINKING NUMBER
#BETWEEN THE LINE p_0-p_1 and the lines p_i-p_i+1 FOR i from 2 to n
#---------------------------------------------------------------------

def seqlkcalc(point_list):

    link_one_rest = []
    if len(point_list)>=4:
        for i in range(2, len(point_list)-1):
            first_link = Segquad([[point_list[0],point_list[1]],[point_list[i],point_list[i+1]]])
            link_one_rest.append(round(first_link.segatan(),3))
    else:
        link_one_rest.append('not enough points!')

    return link_one_rest

#---------------------------------------------------------------------
#FIND AND PARSE ALL CIF FILES IN FOLDER.
#ITERATE THROUGH EACH FILE AND PRODUCE:
# - FRACTIONAL CO-ORDINATES FOR ATOMS
# - CENTROID CO-ORDINATES FOR EACH MOLECULE IN A UNIT CELL
# - SOME LIST OF LINKING NUMBERS
#---------------------------------------------------------------------
files = list(gemmi.CifWalk(inputdir))
print(files)

#Creates list of file names

file_list = []
link_list = []

for file_path in files:
    file_path_split = file_path.split("/")
    file_name = file_path_split[len(file_path_split)-1]
    file_list.append(file_name)
    doc = cif.read(file_path)
    block = doc[0]

    for b in doc:
        if (b.find_loop('_atom_site_') != None):
            block = b

    print("Processed file --> " + file_name)
    crystal_reader = CrystalReader(file_path)
    crystal = crystal_reader[0]
    crystal.assign_bonds()
    packed_molecules = crystal.packing(box_dimensions = ((0,0,0),(1,1,1)), inclusion='CentroidIncluded')
    packed_molecules.normalise_labels()
  
    
    adta_molecules = []
    cent = []
    cent_points = []
    for comp in packed_molecules.components:
        if (len(comp.atoms) > 1):
            adta_molecules.append(comp)
            cent.append(MD.atom_centroid(*list(a for a in comp.atoms)))


    for c in cent:
        cent_points.append([round(c[0],3), round(c[1], 3),round(c[2],3)])

    link_list.append(seqlkcalc(cent_points))

print(link_list)

link_dic = dict(zip(file_list,link_list))

with open(outputdir + '/lkno.csv', 'w') as f:
    for key in link_dic.keys():
        f.write("%s,%s\n"%(key,link_dic[key]))
        
#---------------------------------------------------------------------
#GIVE PAIRWISE RMSD AND POWDER COMPARISONS BETWEEN ALL CIF FILES
#IN INPUT DIRECTORY. ALSO GIVE ABSOLUTE DIFFERENCE OF THE TWO MAXIMUM
#LINKING NUMBERS CALCULATED FOR EACH FILE BY BY THE SEQLINKCALC
#FUNCTION 
#---------------------------------------------------------------------
    

full_list = []
for first_comp in files:
    i = files.index(first_comp)
    rest_list = files[i+1:]
    for second_comp in rest_list:
            first_comp_split = first_comp.split("/")
            second_comp_split = second_comp.split("/")
            first_file_name = first_comp_split[len(first_comp_split)-1]
            second_file_name = second_comp_split[len(second_comp_split)-1]
            one_crys_read = CrystalReader(first_comp)
            two_crys_read = CrystalReader(second_comp)
            one_crys = one_crys_read[0]
            two_crys = two_crys_read[0]
    
            powder_sim = CD.PowderPattern.from_crystal(one_crys)
            powder_comp = CD.PowderPattern.from_crystal(two_crys)

            comp = similarity_engine.compare(one_crys,two_crys)
            powd = powder_sim.similarity(powder_comp)

            try:
                lkcomp = round(abs(max(link_dic[first_file_name]) - max(link_dic[second_file_name])),3)
            except TypeError:
                lkcomp = 'n/a'

            if powd is None:
                pow_input = ' n/a '
            else:
                pow_input = round(powd,3)

            if comp is None:
                comp_input = 'n/a'
            else:
                comp_input = round(comp.rmsd, 3)

            full_list.append([first_file_name, second_file_name, comp_input, pow_input, lkcomp])

            
            with open(outputdir + '/powdercomp.csv','w',newline = '') as p:
                powdwriter = csv.writer(p,delimiter=',',
                            quoting=csv.QUOTE_NONNUMERIC)
                powdwriter.writerows(full_list)

