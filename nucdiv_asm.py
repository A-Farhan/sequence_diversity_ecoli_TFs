### Calculate nucleotide diversity from MSA of refseq(assembled) data

# require
import os, sys
from utilities3.basic import *
from utilities3.custom_classes import *
from joblib import Parallel, delayed

# start clock
start = time()

# number of threads
nt = int(sys.argv[1])

### Generate instances of class MSA ###
# paths
dr = 'codon_msa'
f_out = 'nucdiv_asm.tab'

# list of cds FASTA files
fls = [ os.path.join( dr, i) for i in os.listdir(dr)]

# Function to tabulate analysis results of multiple MSAs
def tabulate_msa_data(fname):
    # check if file exists and is non-zero
    if not filetest_s(fname): return
    msa = MSA(fname,seqtype='n')
    out = [ msa.gene, msa.ncol, msa.nogapcol, msa.ndif, msa.div()]
    return out

# header for the table
header = [ 'gene', 'nuc_len', 'nuc_nogap', 'nuc_ndif', 'nuc_div']

# execute, if output doesn't exist already
if not os.path.exists(f_out):
    print("%s: generating genes' summaries on nucleotide sequence\n"%timer(start))
    stat = filter( None, Parallel(n_jobs=nt)( delayed(tabulate_msa_data) (i) for i in fls) )
    stat = { i[0]:i[1:] for i in stat}
    print("%s: Assembling data into a table and writing to file\n"%timer(start))
    divtab = [ [k] + stat[k] for k in sorted(stat.keys())]
    writecsv( data=divtab, fl=f_out, sep='\t',header=header)
    print("%s: program finished\n"%timer(start))
else:
    print("%s: output file %s already exists\n"%(timer(start),f_out)) 

