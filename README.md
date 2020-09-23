# Powder_RMSD_Link_Compare
Make RMSD, xPRD and if required 'Linking number' comparisons between a selection of CIF files in a folder

Current default setting produces only RMSD and xPRD comparison tables

Relies on CCDC package and must therefore be run from within the CSD Python API.
Requires installation of pandas,sympy and gemmi

Runs with FOUR arguments giving an input and output directory (in that order!) plus a start and end index for the files in the input directory

Input directory must contain one or more CIF files

Output directory will contain comparison files rmsd_compare.csv and powder_compare.csv

In both cases there are 'zero' entries for symmetrical comparisons. 

The current versions of the csv files give comparisons between the first 100 crystals of the T2 dataset.

***HEADER NAMES FOR THE CSV FILES ARE HARDCODED TO WORK WITH T2 DATASET FILENAMES. 
   TO STRIP OUT ONLY THE FILE EXTENSION  LINE 207 SHOULD READ
   
   filename_list = [str(i[len(i)-1][:-4]) for i in filename_split]
