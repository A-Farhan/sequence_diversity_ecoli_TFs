# for a given set of sequencing runs, generate a table of variants with their frequencies
# excluding ancestral

# require
import os, sys, re
from collections import Counter
from utilities3.basic import *
from utilities3.specific import *
from utilities3.lists_dicts import *
from utilities3.Data import seqrecord 

# arguments
wdr = sys.argv[1] # full path

fv = 'variants.tab'
fout = os.path.join( wdr, 'vars_freqs.tab')

# convert seqrecord to a dict of codon lists
seqrecord = {k:n2c(str(v))[:-1] for k,v in seqrecord.items()}

# codon table
ct = codon_table[11].forward_table
stops = codon_table[11].stop_codons

# function to get a list of non-syn variants given a group of runs
def ls_nsv_grp(dr):
    # list of sequencing runs
    runs = [ i for i in os.listdir(dr) if 'RR' in i]
    # initialize list to hold variants
    allv = []
    # for every run
    for seqrun in runs:
        # input file of variants
        fin = os.path.join( dr, seqrun, fv)
        # list of variants
        variants = [ i for i in readcsv(fin,sep='\t')[1:]]
        # initialize list to hold genes with nonsense mutations
        nsg = []
        # initialize list to hold selected variants
        sel = []
        # for every variant
        for v in variants:
            # if the variant is in a gene with a nonsense var,skip
            if v[0] in nsg: continue
            # skip variant, if non-sense
            if v[2] in stops: 
                nsg.append(v[0])
                continue
            # skip variant, if gapped
            if v[2] == '.': continue
            # 0-based codon position
            p = int(v[1])-1
            # skip variant, if position exceeds ref
            if p >= len(seqrecord[v[0]]): continue
            # reference codon
            r = seqrecord[v[0]][p]
            # skip if var is synonymous
            if ct[v[2]] == ct[r]: continue
            # else
            sel.append(''.join(v))
        # add to the pool
        allv.extend(sel)
    # sort variants in descending order of frequency
    out = [ [ i[0], i[1]] for i in Counter(allv).most_common()]
    # remove variants present in all runs ( ancestral)
    out = [ i for i in out if i[1] != len(runs)]
    return out

## main
if not os.path.exists(fout):
    vftab = ls_nsv_grp(wdr)
    writecsv( data=vftab, fl=fout, sep='\t', \
        header=['variants','frequency'])
