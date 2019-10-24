# select a set of sequencing runs satisfying minimum distance threshold given a distance matrix

# require
import os, sys
import numpy as np
from utilities3.basic import *

# start clock
start = time()

# paths
pdr = sys.argv[1]

# arguments
cut = float(sys.argv[2]) # minimum distance threshold

# files
fdm = os.path.join( pdr, 'dmat_cdis.tab')
fout = os.path.join( pdr, 'selected_runs.txt')

# load distance matrix
dmat = readcsv(fdm,sep='\t')
# list of sequencing runs
names = dmat[0]
# array form of distance matrix
dmat = np.array(dmat[1:],dtype='float')

# number of runs
L = len(names)

if not os.path.exists(fout):
    ## execute
    sel = selbydist( dmat, cut, names)
    print("%s: set of runs are now selected\n"%timer(start))

    # write to file
    with open(fout,'w') as flob:
        for i in sel: flob.write( i + '\n')
else:
    print("%s: selection was made already\n"%timer(start))

