## Make a table of genes with hits passing a threshold of gaps

# require
import os, sys
from utilities3.basic import *

# maximum gap
maxg = int(sys.argv[1])

# output file
fout = 'gene_ortho.tab'

# list of table files
fls = [ i for i in os.listdir() if i.endswith('_aln.tab')]

out = []
for f in fls:
    # gene name
    gene = f.replace('_aln.tab','')
    # load data
    data = readcsv(f,sep='\t')[1:]
    # subset data by gap threshold
    sub = [ i for i in data if float(i[-1]) <= maxg]
    # number of hits
    nh = len(data)
    # number of selected hits
    nsh = len(sub)
    out.append([gene,nh,nsh])
# write output table to file
writecsv(data=out,fl=fout,sep='\t',header=['gene','total','selected'])
