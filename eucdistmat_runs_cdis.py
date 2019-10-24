# generate a euclidean distance matrix from codon distances of runs from reference on 125 genes

# require
import os, sys
import numpy as np
from scipy.spatial import distance
from utilities3.basic import *

# paths
# working directory
wdr = sys.argv[1]
os.chdir(wdr)

# files
fin = 'codon_distance.tab'
fout = 'dmat_cdis.tab'

# load data extract runs & genes
data = readcsv( fin, sep='\t')
runs = data[0][1:]
data = data[1:]
tdata = list(zip(*data))
genes = tdata[0]
cdmat = np.array( tdata[1:], dtype='float')

# no. of runs and genes
R,G = cdmat.shape

# initialize a distance matrix
dmat = np.array( [ [0 for c in range(R)] for r in range(R)], \
        dtype='float')
# compute euclidean distances
for x in range(R-1):
    for y in range(x+1,R):
        d = distance.euclidean( cdmat[x,], cdmat[y,])
        dmat[x,y] = dmat[y,x] = round(d,2)
# write output to file
writecsv(data=dmat, fl=fout, sep='\t', header = runs)
