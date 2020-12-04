#---------------------------------------------------------------------
#PACKAGE IMPORTS
#---------------------------------------------------------------------

import sys
import os
import ccdc
import random
import time
import matplotlib.pyplot as plt
from ccdc import search, io, molecule
from ccdc.io import MoleculeReader, CrystalReader, EntryReader, CrystalWriter
from ccdc.descriptors import MolecularDescriptors as MD
from ccdc.descriptors import CrystalDescriptors as CD
from ccdc.crystal import PackingSimilarity
similarity_engine = PackingSimilarity()
import numpy as np
import gemmi
from gemmi import cif
import csv
import itertools

#---------------------------------------------------------------------
#SYSTEM PARAMETERS
#---------------------------------------------------------------------

# Check parameters
if (len(sys.argv) != 3):
    print("Must have 2 parameters: input dir!")
    sys.exit()
inputdir = sys.argv[1]
outputdir = sys.argv[2]

# Check they are directories
if not(os.path.isdir(inputdir) & os.path.isdir(outputdir)) :
    print("Parameters are not directories!")
    sys.exit()
        
#---------------------------------------------------------------------
#EXTRACT n RANDOM FILES FROM THE SOURCE
#---------------------------------------------------------------------
def file_select(n):

    files = list(gemmi.CifWalk(inputdir))

    file_select = random.sample(range(0, len(files)-1), n)

    sorted_index_list = [files[i] for i in file_select]
                
    print("Files to be compared:")
    print(sorted_index_list)

    return sorted_index_list

#---------------------------------------------------------------------
#EXTRACT K RANDOM FILES AND DO RMSD COMPARISONS M TIMES
#OUTPUT PLOT OF AVERAGE TIMES TAKEN AND TOTAL AVERAGES
#---------------------------------------------------------------------

def rmsdtime(k,m):

    rpt = 0
    xax = []

    process_time_list = []
    compare_time_list = []
    total_time_list = []

    while rpt <=m:



        file_extract = file_select(k)

        filename_split = [i.split("/") for i in file_extract]
        filename_list = [str(i[len(i)-1][4:-4]) for i in filename_split]
        print("Files to be compared:")
        print(filename_list)


        a = time.perf_counter()
        crystals = [CrystalReader(c)[0] for c in file_extract]
        #packed_crystals = [c.packing() for c in crystals]
        b = time.perf_counter()

        process_time_list.append(b-a)

        print(f'Processing time (s) for %s molecules:', k)
        print((b-a))
        print(process_time_list)
    

        c = time.perf_counter()
        for i, j in itertools.combinations_with_replacement(range(len(crystals)), 2):

            comp = similarity_engine.compare(crystals[i],crystals[j])
        d = time.perf_counter()

        print(f'Comparison time (s) for %s molecules:', k)
        print((d-c) + (a-b))

        compare_time_list.append(d-c)
        print([compare_time_list])

        total_time_list.append((b-a) + (d-c))
        print([total_time_list])

        rpt+=1
        xax.append(rpt)
        

    plt.plot(xax, process_time_list, label = 'Processing Time')
    plt.plot(xax, compare_time_list, label = 'Comparison Time')
    plt.plot(xax, compare_time_list, label = 'Total Time')
    plt.xlabel('Run number')
    plt.ylabel('Computation time')
    plt.title(f'RMSD Comparison')
    plt.legend()

    plt.savefig(outputdir+ "/RMSD_time_test.png")
    plt.cla()    

    print(f'Average rmsd processing time (s) for 5 molecules:')
    print(sum(process_time_list)/len(process_time_list))
    print(f'Average rmsed pairwise comparison time for 5 molecules')
    print(sum(compare_time_list)/len(compare_time_list))
    print(f'Total rmsd time for %s molecules', k)
    print(sum(total_time_list)/len(total_time_list))

    return [sum(process_time_list)/len(process_time_list), sum(compare_time_list)/len(compare_time_list), sum(total_time_list)/len(total_time_list)]

#---------------------------------------------------------------------
#EXTRACT K RANDOM FILES AND DO POWDER COMPARISONS M TIMES.
#OUTPUT AVERAGE TIME PLOT AND OVERALL AVERAGE
#---------------------------------------------------------------------
def powdertime(k,m):

    rpt = 0
    xax = []

    process_time_list = []
    compare_time_list = []
    total_time_list = []

    while rpt <=m:



        file_extract = file_select(k)

        #filename_split = [i.split("/") for i in file_extract]
        #filename_list = [str(i[len(i)-1][4:-4]) for i in filename_split]

        a = time.perf_counter()
        crystals = [CrystalReader(c)[0] for c in file_extract]
        #packed_crystals = [c.packing() for c in crystals]
        powders  = [CD.PowderPattern.from_crystal(c) for c in crystals]
        b = time.perf_counter()

        process_time_list.append(b-a)

        print(f'Processing time (s) for 5 molecules')
        print((b- a))
    

        c = time.perf_counter()
        for i, j in itertools.combinations_with_replacement(range(len(crystals)), 2):
            powd = powders[i].similarity(powders[j])
        d = time.perf_counter()

        compare_time_list.append(d-c)
        total_time_list.append((b-a) + (d-c))

        print(f'Compariosn time (s) for 5 molecules:')

        rpt+=1
        xax.append(rpt)
        
        
        
    plt.plot(xax, process_time_list, label = 'Processing Time')
    plt.plot(xax, compare_time_list, label = 'Comparison Time')
    plt.plot(xax, total_time_list, label = 'Total Time')
    plt.xlabel('Run number')
    plt.ylabel('Computation time')
    plt.title('Powder Comparison')
    plt.legend()

    plt.savefig(outputdir+ "/powder_time_test.png")
    plt.cla()  

    print(f'Average powder processing time (s) for %s molecules:', k)
    print(sum(process_time_list)/len(process_time_list))
    print(f'Average powder pairwise comparison time for %s molecules', k)
    print(sum(compare_time_list)/len(compare_time_list))
    print(f'Total powder time for %s molecules', k)
    print(sum(total_time_list)/len(total_time_list))

    return [sum(process_time_list)/len(process_time_list), sum(compare_time_list)/len(compare_time_list), sum(total_time_list)/len(total_time_list)]



