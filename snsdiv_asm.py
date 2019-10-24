### Calculate synonymous and non-synonymous diversity from MSA of refseq(assembled) data

# require
import os, sys
from utilities3.basic import *
from utilities3.specific import *
from utilities3.custom_classes import *
from utilities3.Data import genes_wo_additionalsites as genes
from joblib import Parallel, delayed
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.special import binom

# start clock
start = time()

# number of threads
nt = int(sys.argv[1])
home = os.path.expanduser('~')
pdr = os.path.join( home, 'varspots')
wdr = sys.argv[2] # path relative to project directory
wdr = os.path.join( pdr, wdr)
os.chdir(wdr)

### Generate instances of class MSA ###
# paths
dr_n = 'codon_msa'
f_out = 'snsdiv_asm.tab'

# list of cds FASTA files
fls_n = [ os.path.join( dr_n, i+'.fasta') for i in genes]

### Auxilliary functions ###

## function to get a dict of genome ids with syn & nonsyn site count
# and the codon sequence
def msasiteseq(msa):
    # list of seq records
    aln = list(msa.aln)
    # array form of MSA without gaps
    arr = msa.remove_gapcols()
    # dimensions
    R,C = arr.shape
    # output dict
    out = {}
    for x in range(R):
        sequence = ''.join(arr[x,])
        Id = aln[x].id
        rec = SeqRecord( Seq(sequence,IUPAC.unambiguous_dna), id=Id)
        s = recsitecount(rec)
        if s is not None:
            out[s[0]] = s[3:5] + [n2c(sequence)]
    return out

## function to get synonymous & non-synonymous diversity 
# using above dict
def snsdiv(siteseqdict):
    dct = siteseqdict
    # list of genomes
    gens = list(dct.keys())
    # no. of genomes
    L = len(gens)
    # no. of unique pairs
    P = binom(L,2)
    # initialize Dn & Ds as zeros
    Dn,Ds = 0,0
    for x in range(0,L-1):
        seq1 = dct[gens[x]][2]
        for y in range(x+1,L):
            S = ( dct[gens[x]][0] + dct[gens[y]][0] )/2
            N = ( dct[gens[x]][1] + dct[gens[y]][1] )/2 
            seq2 = dct[gens[y]][2]
            res = DnDs(seq1,seq2)
            Dn += res[0]/N
            Ds += res[1]/S
    sdiv = round( Ds/P, 4)
    nsdiv = round( Dn/P, 4)
    return [sdiv,nsdiv]

## function to tabulate analysis results of multiple MSAs
def tabulate_msa_data(fname,seqtype):
    # check if file exists and is non-zero
    if not filetest_s(fname): return
    msa = MSA(fname,seqtype)
    gene = msa.gene
    print("%s: For %s\n"%(timer(start),gene))
    siteseq = msasiteseq(msa)
    R = len(siteseq)
    if R <= 1: return

    Sd,Nd = snsdiv(siteseq)
    out = [gene,R,Sd,Nd]
    return out

# header for the table
header = [ 'gene', 'runs' ,'syn', 'nonsyn']

# execute, if output doesn't exist already
if not os.path.exists(f_out):
    print("%s: generating genes' summaries on nucleotide sequence\n"%timer(start))
    out = filter(None,Parallel(n_jobs=nt)( delayed(tabulate_msa_data) \
            (i,seqtype='n') for i in fls_n ))
    writecsv( data=out, fl=f_out, sep='\t',header=header)
    print("%s: program finished\n"%timer(start))
else:
    print("%s: output file %s already exists\n"%(timer(start),f_out)) 
