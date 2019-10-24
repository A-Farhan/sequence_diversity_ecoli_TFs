# make a table of codon distances from sequencing runs in a given directory
# exclude runs without input files or not enough genes

# require
import os, sys
from utilities3.basic import * 
from utilities3.wgs_diversity import rungenedict
from utilities3.Data import rawhdgenes as genes

# working directory
wdr = sys.argv[1]
os.chdir(wdr)

# file name for list of codon distances
fcd = 'codon_distances.txt'
# output file path
fout = 'codon_distance.tab'

# minimum gene threshold
minG = 3000

# list of runs which have input file
runs = [ i for i in os.listdir() if 'RR' in i and \
        filetest_s( os.path.join( i, fcd))]

# dict of runs with their genes
run_genes = rungenedict(runs)

# exclude runs without sufficient genes
runs = [ i for i in runs if len( run_genes[i]) > minG]

# table of distances
tab = []
for r in runs:
    fin = os.path.join( r, fcd)
    if not os.path.exists(fin): continue
    with open(fin) as flob:
        data = [ float(i.strip('\n')) for i in flob]
        tab.append([r]+data)

tab = [ ['runs'] + genes ] + tab 
# transpose the above table
ttab = [ i for i in zip(*tab)]
# write to file
writecsv( data=ttab, fl=fout, sep='\t')
