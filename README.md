# Powder_RMSD_Link_Compare
Make RMSD, xPRD and 'Linking number' comparisons between all CIF files in a folder

Relies on CCDC package and must therefore be run from within the CSD Python API.
Requires installation of sympy and gemmi

Runs with two arguments giving an output and input directory (in that order!)

Output directory must contain one or more CIF files

Input directory will contain the following files:

lkno.csv gives a table of all CIF filenames in the output folder along with, currently, a vector giving the linking number between a line       joining the first two centroids in a crystal packing and a line between centroids 2&3, 3&4 etc. If there are only three centroids in the crystal it outputs 'not enough points!

powdercomp.csv gives a pairwise comparison between each file in the folder. The first two entries of each row are the files being compared, the second the RMSD, the third the powder comparison. The fourth is simply the absolute difference between the maximum linking number in each link

The current versions of the csv files give comparison between experimental ADTA structures taken from Cui et al. 'Mining predicted crystal structure landscapes with high throughput crystallisation: old molecules, new insights' (DOI: 10.1039/c9sc02832c)
