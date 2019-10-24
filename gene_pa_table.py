# make a table of gene presence/absence in different runs of a dataset
## the table includes all runs with gene lists irrespective of the number of genes in them
## only the runs without genelists are excluded

# require
import os, sys
from utilities3.basic import *
from utilities3.specific import *

# start clock
start = time()

# full path to the data directory
wdr = sys.argv[1]

# files
f_pa = 'genes_present.txt'
f_out = os.path.join( wdr, 'genes_pa.tab')

# list of runs with gene lists
# excluding those that don't have the input file
runs = [ i for i in os.listdir(wdr) if filetest_s(os.path.join(wdr,i,f_pa))]

# initialize a list to hold all genes present in the dataset
genes = []
# go through each run directory
for r in runs:    
    fl = os.path.join( wdr, r, f_pa)
    # load the list of genes inferred present
    with open(fl) as flob:
        data = [ i.strip('\n') for i in flob]
    # go through each gene
    for gene in data:
        # if the gene is not already present in the above list of genes
        if gene not in genes:
            # then include it
            genes.extend([gene])
print( "%s: list of genes present in the dataset created\n"%timer(start))

# no. of runs and genes
R = len(runs)
G = len(genes)

# initialize the output table
tab = [ [ genes[y]] + [ 0 for x in range(R)] \
        for y in range(G)]

# go through each run directory again
for x, r in enumerate(runs):
    # load the run's gene list
    fl = os.path.join( wdr, r, f_pa)
    with open(fl) as flob:
        data = [ i.strip('\n') for i in flob]
    
    # go through the list of all genes
    for y in range(G):
        # if the gene is present in the data
        if genes[y] in data:
            tab[y][x+1] = 1
print( "%s: output table generated\n"%timer(start))

# order the output table
## st more filled rows are above
tab = order_pamat(tab)
print( "%s: output table ordered\n"%timer(start))

# write the table to file
writecsv( data = tab, fl = f_out, sep = '\t', \
            header = ['genes']+runs)
print( "%s: table saved as %s\n"%(timer(start),f_out))
