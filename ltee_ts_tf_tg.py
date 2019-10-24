# process AF time series data from ltee classifying genes as tf & targets
# keep generation info 

## the script doesn't check if output file already exists

# require
import os, re, sys
from utilities3.basic import *
from utilities3.Data import get_trn
from utilities3.lists_dicts import flatten

# pass full path to data directory and population code as command-line argument
dr = sys.argv[1]
pop=sys.argv[2]

# paths
fin = os.path.join( dr, pop+'_annotated_timecourse.txt')
fout = os.path.join( dr, pop+'_ts.tab')

# TRN
trn = get_trn()
tfs = list(trn.keys())
targets = list( set( i for i in flatten(trn.values()) ))

# load time series data
data = readcsv(fin,sep=',')
header = data[0]
data = data[1:]

# keep passed and missense mutations
data = [ [ j.strip(' ') for j in i] for i in data \
        if i[12] == ' PASS' and i[3] == ' missense']

# indices of columns of generations 
# starting with AC
gencols = [ i for i in header if bool(re.search('AC',i))]
gencols = [ i for i in gencols if int(i.split(':')[1])<=60000]
ixs = [ header.index(i) for i in gencols]

# get generation column names
gens = [ re.sub(pattern=' AC:',repl='',string=i) for i in gencols] 

# extract data on mutated gene, annotation and frequencies
fdata = [ i[:2] + [ round(float(i[x])/float(i[x+1]),2) \
        if i[x+1] != '0.0' else 0 for x in ixs] for i in data]

# identify genes as tfs or targets and keep only these mutations
out = []
for row in fdata:
    if row[1] in tfs:
        bit = 1
    elif row[1] in targets:
        bit = 0
    else: continue
    entry = [ bit, row[1], row[0]] + row[2:]
    out.append(entry)

# prepare output header
H = [ 'class', 'gene', 'position'] + gens 
# write table to a new file
writecsv(data=out,fl=fout,sep='\t',header=H)
