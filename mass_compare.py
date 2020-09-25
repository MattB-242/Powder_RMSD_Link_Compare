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
import itertools
import seaborn as sns

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
        
#---------------------------------------------------------------------
#EXTRACT LIST OF FILES SORTED BY LOWEST ENERGY AND MAKE LIST
#OF FILES IN POSITIONS st TO e
#---------------------------------------------------------------------

files = list(gemmi.CifWalk(inputdir))

energy_index = pd.read_csv(inputdir + '/E_D_List.csv', header = 0)
energy_index.sort_values(by = ['ENERGY'])

sorted_full_list = energy_index['ID'].tolist()
sorted_index_list = sorted_full_list[st:e]

print(sorted_index_list)

file_extract = [inputdir + "/" + g+'.cif' for g in sorted_index_list if [g in f for f in files]]
                
print (file_extract)

#---------------------------------------------------------------------
#FOR FILES IN THE EXTRACTED LIST
#EXPORT PAIRWISE COMPARISON TABLES FOR RMSD, RMSD
#IF ALL SET TO TRUE.

#DEFAULT SETTINGS RUN POWDER AND RMSD FOR ALL FILES IN DIRECTORY
#---------------------------------------------------------------------

def molcomp(start, end, pcomp=True, rmsd_threshold=2, rounding=6, rcomp=True, molcomp = True):

    filename_split = [i.split("/") for i in file_extract]
    filename_list = [str(i[len(i)-1][4:-4]) for i in filename_split]

    l = len(file_extract)
    if molcomp:
        molarray = np.zeros((l, l))

    if pcomp:
        powdarray = np.zeros((l, l))

    if rcomp:
        rmsdarray = np.zeros((l, l))

    crystals = [CrystalReader(c)[0] for c in file_extract]
    powders  = [CD.PowderPattern.from_crystal(c) for c in crystals]

    for i, j in itertools.combinations_with_replacement(range(len(crystals)), 2):

        comp = similarity_engine.compare(crystals[i],crystals[j])
        powd = powders[i].similarity(powders[j])

        if molcomp:
            try:
                molarray[j][i] += (comp.nmatched_molecules)
                molarray[i][j] += (comp.nmatched_molecules)
            except TypeError:
                molarray[j][i] += 99
                molarray[i][j] += 99

        
        if pcomp:
            if powd is None:
                powdarray[j][i] += 99
                powfarray[i][j] += 99
            else:
                powdarray[j][i] += (round(powd,6))
                powdarray[i][j] += (round(powd,6))

        if rcomp:
            if comp is None:
                rmsdarray[j][i]+= 99
                rmsdarray[i][j]+=99
            elif comp.nmatched_molecules < rmsd_threshold:
                rmsdarray[j][i]+= 88
                rmsdarray[i][j]+= 88
            else:
                rmsdarray[j][i]+=(round(comp.rmsd, 6))
                rmsdarray[i][j]+=(round(comp.rmsd, 6))
            
    if molcomp:
        molframe = pd.DataFrame(molarray, index=filename_list, columns = filename_list)
        molframe.replace(99, 'nan')
        molframe.to_csv(outputdir+"/rmsd_molecules_matched.csv", index = True, header = True, sep = ',')

    if pcomp:
        powdframe = pd.DataFrame(powdarray, index=filename_list, columns = filename_list)
        powdframe.to_csv(outputdir+"/powder_comparison.csv", index = True, header = True, sep = ',')

    if rcomp:
        rmsdframe = pd.DataFrame(rmsdarray, index=filename_list, columns = filename_list) 
        rmsdframe.to_csv(outputdir+"/rmsd_comparison.csv", index = True, header = True, sep = ',')


#-----------------------------------------------------------------------
#RUN LKCOMP OUTPUTS WITH DEFAULT OUTPUT SETTINGS ON FIRST FILES st to e
#-----------------------------------------------------------------------
molcomp(st,e)
