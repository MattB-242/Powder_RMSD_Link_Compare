# Powder_RMSD_Link_Compare
Make RMSD, pXRD comparisons between a selection of CIF files in a folder

Current default setting produces only RMSD and xPRD comparison tables

Relies on CCDC package and must therefore be run from within the CSD Python API.

Requires installation of pandas,sympy and gemmi

Runs with FOUR arguments:

linkcomp.py <input directory> <output directory> <start file index> <end file index>

Input directory must contain one or more CIF files

Output directory will contain comparison files rmsd_compare.csv and powder_compare.csv after run

In both cases there are 'zero' entries for symmetrical and self-comparisons. 

Zero entries in the RMSD comparison may indicate a comparison failure (insufficient molecules matched)

The current versions of the csv files give comparisons between the first 20 crystals of the T2 dataset.

***HEADER NAMES FOR THE CSV FILES ARE HARDCODED TO WORK WITH T2 DATASET FILENAMES. 
   TO STRIP OUT ONLY THE FILE EXTENSION  LINE 98 SHOULD READ
   
   filename_list = [str(i[len(i)-1][:-4]) for i in filename_split]
