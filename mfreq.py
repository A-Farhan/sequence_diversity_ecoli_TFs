## Get frequencies with which TFs and targets are mutated 

## require
import os, sys, re
from utilities3.basic import *
from utilities3.Data import get_trn, genesitecount
from utilities3.lists_dicts import flatten
from collections import Counter

## paths
wdr = sys.argv[1]
fv = os.path.join( wdr, 'vars_freqs.tab')
fout = os.path.join( wdr, 'mfreq.tab')

# TRN
trn = get_trn()
tfs = list(trn.keys())
targets = list( set( i for i in flatten(trn.values()) ))

# nonsyn site count for tfs and targets
sites = { k:v[2] for k,v in genesitecount().items()}
sf = round( sum( [ v for k,v in sites.items() if k in tfs]))
sg = round( sum( [ v for k,v in sites.items() if k in targets]))

# load variant frequency table
variants = [ i[0] for i in readcsv(fv,sep='\t')[1:]]

# extract gene names from variants
mgenes = [ re.split('\d+',i)[0] for i in variants]

# get no. of mutations in each gene
gene_mut = { k:v for k,v in Counter(mgenes).items()}

# mutation per site for TFs and targets
mf = sum([v for k,v in gene_mut.items() if k in tfs])
mg = sum([v for k,v in gene_mut.items() if k in targets])

## make contingency table
ctab = [ [ 'tf', mf, sf - mf], [ 'tg', mg, sg - mg]]
writecsv(data=ctab,fl=fout,sep='\t',\
        header = ['class', 'mutated','not-mutated'])

## ODDs ratio
mpsf = mf/(sf-mf)
mpsg = mg/(sg-mg)
R = mpsf/mpsg
print("Odds Ratio for observed mutation frequency in TFs = %f\n"%R)
