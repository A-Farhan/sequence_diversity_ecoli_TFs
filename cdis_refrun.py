# script to generate codon-distances from reference for a list of runs
# on a list of selected genes 
# sequencing runs which don't have the variants table or don't have sufficient number of genes are excluded

# require
import os, sys
import numpy as np
from joblib import Parallel,delayed
from utilities3.basic import * 
from utilities3.specific import c2n
from utilities3.wgs_diversity \
        import rungenedict,genevardict
from utilities3.Data import seqrecord, rawhdgenes as genes

# start clock
start = time()

# arguments
# number of threads
nt = int(sys.argv[1])

# working directory
wdr = sys.argv[2]
os.chdir(wdr)

# file name for table of variants
fv = 'variants.tab'
# file name for list of codon distances
fcd = 'codon_distances.txt'

# minimum gene threshold
minG = 3000

# list of runs
# excluding those lacking variant file 
runs = [ i for i in os.listdir() if 'RR' in i and filetest_s(os.path.join(i,fv))]

# dict of runs with their genes
run_genes = rungenedict(runs)
print("%s: generated dict of runs with their genes\n"%timer(start))

# exclude runs not satisfying minimum gene threshold
runs = [ i for i in runs if len(run_genes[i]) > minG]

# hold data from multiple run's variant files into a dict
# limit to list of selected genes
gene_vars = { k:v for k,v in genevardict(runs).items() if k in genes}
print("%s: generated dict of gene with their variants\n"%timer(start))

# function to calculate codon-distance of a sequencing run 
# from the reference given gene-run-variant dict
def cdis_ref(seqrun,gene_vars,seqrecord):
    # list of genes mutated in the run
    mgenes = [ k for k,v in gene_vars.items() \
            if seqrun in v.keys()]
    # initialize output dict
    out = {}
    # for every mutated gene
    for g in mgenes:
        # list of variant codons
        varc = list( gene_vars[g][seqrun].values())
        # no. of variant codons
        N = len(varc)
        # no. of codons in the gene
        L = len(seqrecord[g])/3
        # gap count
        gc = varc.count('.')
        # codon-distance 
        cdis = round( (N-gc)/(L-gc), 4)
        out[g] = cdis
    return out

## wrapper function to generate a list of pdistances for given run
def wrapper(seqrun,genes=genes,gene_vars=gene_vars,run_genes=run_genes,seqrecord=seqrecord):
    # output file path
    fout = os.path.join( seqrun, fcd)
    # skip if already exists
    if os.path.exists(fout): return
    print("%s: For %s\n"%(timer(start),seqrun))
    # dict of mutated genes with codon distances from the reference
    dct = cdis_ref(seqrun,gene_vars,seqrecord)
    # initialize output list
    out = []
    # for every gene of interest
    for g in genes:
        # if gene is mutated, then take its codon distance
        if g in dct.keys(): v = dct[g]
        # if gene is not mutated but present in the sequencing run, consider codon distance as zero 
        elif g in run_genes[seqrun]: v = 0 
        #  if neither of the above, use filler value -1 
        else: v = -1
        out.append(v)
    # write output to file
    with open(fout,'w') as flob:
        for i in out: flob.write(str(i)+'\n')

## execute
Parallel(n_jobs=nt) ( map( delayed(wrapper), runs))
print(timer(start))
