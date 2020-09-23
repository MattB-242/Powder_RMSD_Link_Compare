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
import pandas as pd
import csv

#---------------------------------------------------------------------
#SYSTEM PARAMETERS
#---------------------------------------------------------------------

# Check parameters
if (len(sys.argv) != 5):
    print("Must have 4 parameters: input dir, output dir, and a start and end index in the file!")
    sys.exit()
inputdir = sys.argv[1]
outputdir = sys.argv[2]
st = int(sys.argv[3])
e = int(sys.argv[4])

# Check they are directories
if not(os.path.isdir(inputdir) & os.path.isdir(outputdir)) :
    print("Parameters are not directories!")
    sys.exit()

print(st+e)

#---------------------------------------------------------------------
#FIND AND PARSE ALL CIF FILES IN FOLDER.
#ITERATE THROUGH EACH FILE AND PRODUCE:
#   - LIST OF ALL FILE NAMES
#   - CENTROID CO-ORDINATES FOR EACH MOLECULE IN A UNIT CELL
#---------------------------------------------------------------------

#Creates list of file paths for all CIF files in input directory
files = list(gemmi.CifWalk(inputdir))


file_list = []

for file_path in files[st:e]:
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
        
#---------------------------------------------------------------------
#FOR FILES start TO end IN THE INPUT DIRECTORY
#EXPORT PAIRWISE COMPARISON TABLES FOR RMSD, POWDER OR BOTH

#DEFAULT SETTINGS RUN POWDER AND RMSD FOR ALL FILES IN DIRECTORY
#---------------------------------------------------------------------
    
def linkcomp(start, end, pcomp=True, rcomp=True):

    #Creates list of file names with 'job_' and '.cif' stripped out
    full_list = []
    file_extract = files[start:end]
    filename_split = [i.split("/") for i in file_extract]
    filename_list = [str(i[len(i)-1][4:-4]) for i in filename_split]

    #Initiate blank arrays for populating with comparisons
    if pcomp:
        powdarray = np.zeros((len(file_extract), len(file_extract)))

    if rcomp:
        rmsdarray = np.zeros((len(file_extract), len(file_extract)))

    #Extract crystal objects to be compared from the two CIF files
    for first_comp in file_extract:
        i = file_extract.index(first_comp)
        rest_list = file_extract[i+1:end]
        new_entry = []
        for second_comp in rest_list:
            j = file_extract.index(second_comp)
            first_comp_split = first_comp.split("/")
            second_comp_split = second_comp.split("/")
            one_crys_read = CrystalReader(first_comp)
            two_crys_read = CrystalReader(second_comp)
            one_crys = one_crys_read[0]
            two_crys = two_crys_read[0]
    
            powder_sim = CD.PowderPattern.from_crystal(one_crys)
            powder_comp = CD.PowderPattern.from_crystal(two_crys)

            #Run comparisons and add to initialised matrices if active
            comp = similarity_engine.compare(one_crys,two_crys)
            powd = powder_sim.similarity(powder_comp)
            
            if pcomp:
                if powd is None:
                    powdarray[j][i] += np.nan
                else:
                    powdarray[j][i]+= (round(powd,3))

            if rcomp:
                if comp is None:
                    rmsdarray[j][i]+= np.nan
                else:
                    rmsdarray[j][i]+=(round(comp.rmsd, 3))

    #Save comparisons as csv with file names as both column and row
    if pcomp:
        powdframe = pd.DataFrame(powdarray, index=filename_list, columns = filename_list)
        powdframe.to_csv(outputdir+"/powder_comparison.csv", index = True, header = True, sep = ',')

    if rcomp:
        powdframe = pd.DataFrame(rmsdarray, index=filename_list, columns = filename_list)
        powdframe.to_csv(outputdir+"/rmsd_comparison.csv", index = True, header = True, sep = ',')


#---------------------------------------------------------------------
#RUN LKCOMP OUTPUTS WITH DEFAULT OUTPUT SETTINGS ON FIRST n FILES
#---------------------------------------------------------------------
linkcomp(st,e)
