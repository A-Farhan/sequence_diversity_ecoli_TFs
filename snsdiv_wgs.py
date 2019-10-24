# Calculate synonymous & non-synonymous diversity from WGS data
# sequencing runs which don't have the variants table or don't have sufficient number of genes are excluded

# require
import sys
from joblib import Parallel, delayed
from utilities3.basic import *
from utilities3.specific import *
from utilities3.wgs_diversity import *
from utilities3.Data import genes_wo_additionalsites as genes

# start clock
start = time()

# number of threads
nt = int(sys.argv[1])
wdr = sys.argv[2] # path 
os.chdir(wdr)

# paths
fv = 'variants.tab'
fdiv = 'snsdiv_wgs.tab'

# minimum gene threshold
minG = 3000

# list of sequencing runs
runs = [ i for i in os.listdir() if 'RR' in i]

# test which runs have the variant file
# exclude those which don't
runs = [ i for i in runs if filetest_s(os.path.join(i,fv))]

# dict of runs with their genes
run_genes = rungenedict(runs)
print("%s: generated dict of runs with their genes\n"%timer(start))

# exclude runs not satisfying minimum gene threshold
runs = [ i for i in runs if len(run_genes[i]) > minG]

# hold data from multiple run's variant files into a dict
gene_vars = genevardict(runs)
print("%s: generated dict of gene with their variants\n"%timer(start))

## wrapper functions
def wrapper(gene,runs=runs,gene_vars=gene_vars,run_genes=run_genes):
    print("%s: For %s\n"%(timer(start),gene))
    # generate a codon matrix from the above gene variant dict
    mat = generunvardict_to_codonmat(gene,runs,gene_vars,run_genes)
    if mat is None: return  
    # compute SNS diversity from the above matrix
    return snsdiv_wgscodonmat(gene,mat)

## apply the wrapper function parallely to all genes with variants
out = Parallel( n_jobs = nt) ( map( delayed(wrapper), genes))
out = filter(None,out)
writecsv(data=out,fl=fdiv,sep='\t',header=['gene','runs','syn','nonsyn'])
print("%s: Output table was generated and written to file\n"%timer(start))
