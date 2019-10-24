## Functions required to calculate nucleotide diversity from WGS data
import os
from utilities3.basic import *
from utilities3.specific import *
from utilities3.Data import seqrecord
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table
import numpy as np
from scipy.special import binom

# Function to hold data from multiple run files into a dict of genes with variants
def genevardict(runlist): # runlist refers to a list of sequencing runs
    # initialize the dict as empty
    out = {}
    # go through each sequencing run
    for seqrun in runlist:
        # file path for variants table of the run
        fpath = os.path.join(seqrun,'variants.tab')
        # if the above file doesn't exist, skip to next run
        if not filetest_s(fpath): continue
        # load data from the above table
        data = readcsv(fpath,sep='\t')[1:]
        # traverse the table row-wise
        for row in data:
            g = row[0] #gene
            # first add the gene as key, if not already present
            if g not in out.keys():
                out[g] = {}
            # then add sequencing run as the second key, if not already present
            if seqrun not in out[g].keys():
                out[g][seqrun] = {}
            # then add zero-indexed position as the third key with the variant codon as value
            out[g][seqrun][int(row[1])-1] = row[2]
    return out

# Function to hold genelist of sequencing runs into a dict
def rungenedict(runlist): # runlist refers to a list of sequencing runs
    out = {}
    for seqrun in runlist:
        fpath = os.path.join(seqrun,'genes_present.txt')
        # skip the run if above file doesn't exist
        if not filetest_s(fpath): continue
        with open(fpath) as flob:
            genelist = [i.strip('\n') for i in flob]
        out[seqrun] = genelist
    return out

# Function to convert a dict of genes with variants across runs to a codon matrix
# removing sites with gaps
# arguments:
# gene: name of the gene, principal mandatory argument
# runs: list of sequencing runs with variant table and gene lists
# gene_vars: dict with genes as 1st key, runs as 2nd, codon position as 3rd and variants as values
# run_genes: dict or runs with their gene lists
# seqrecord: dict of sequence records of reference genes
# stops: list of stop codons

def generunvardict_to_codonmat(gene,runs,gene_vars,run_genes,seqrecord=seqrecord,stops=codon_table[11].stop_codons):
    # no. of runs
    N = len(runs)
    # nucleotide sequence
    nseq = str(seqrecord[gene])
    # length of nucleotide sequence
    ln = len(nseq)
    # codon sequence
    cseq = n2c(nseq)
    # length of codon sequence
    lc = len(cseq)
    # last codon position's index
    spos = lc-1
    
    # initialize codon matrix
    mat = np.array( [ cseq for r in range(N)])
    # dimensions of the above matrix
    R,C = mat.shape
 
    # if the gene is mutated in any of the sequencing runs
    if gene in gene_vars.keys():
        # dict of variants
        vardict = gene_vars[gene]
        # mutated positions 
        pos = sorted( list( set( k for v in vardict.values() for k in v.keys())))
        # terminate the function if the reported mutated position is not compatible with gene length
        if pos[-1] >= spos: return
    
        # identify the runs with mutations at the stop position to non-stop codons
        stopruns = [ k for k,v in vardict.items() if spos in v.keys() and v[spos] not in stops]

        # replace ref codons with variants at mutated positions
        for r,seqrun in enumerate(runs):
            if seqrun in vardict.keys():
                for p in pos:
                    if p in vardict[seqrun].keys():
                        mat[r,p] = vardict[seqrun][p]   
    
        # remove runs with following issues
        rmrows = [ r for r in range(R) if \
                    # gene absent from the run
                    gene not in run_genes[runs[r]] \
                    # intermediate stop codons
                    or any( i in stops for i in mat[r,:-1]) \
                    # mutation in the stop codon
                    or runs[r] in stopruns]

    # if gene is NOT mutated 
    else:
        # indices of runs which don't have the gene
        rmrows = [ r for r in range(R) if gene not in run_genes[runs[r]] ]
    
    # remove runs identified by above indices
    mat = np.delete(mat,rmrows,0)
    
    # if no sequencing run is left, return None
    if mat.shape[0] == 0: return
           
    # remove columns with gaps
    rmcols = [ c for c in range(C) if '.' in mat[:,c]]
    mat = np.delete(mat,rmcols,1)

    # remove stop codon
    mat = mat[:,:-1]
    return mat

# Function to compute nucleotide divesity given a codon matrix
def nucdiv_wgscodonmat(gene,mat):
    # number of sites are 3 times no. of columns ( codon to nucleotide)
    nsites = mat.shape[1]*3
    # convert codon matrix to nucleotide matrix
    nmat = np.apply_along_axis(c2n,1,mat)
    # indices of invariant columns
    rmcols = [ c for c in range(nmat.shape[1]) if len(set(nmat[:,c]))==1]
    # if all columns are invariant
    if len(rmcols) == nmat.shape[1]:
        out = [gene, nmat.shape[0], 0, 0]
    else:
        nmat = np.delete(nmat,rmcols,1)
        R,C = nmat.shape
        S = sum( sum( nmat[x,] != nmat[y,]) for x in range(0,R-1) for y in range(x+1,R))
        D = round( S/(nsites*binom(R,2)), 4)
        out = [gene,R,C,D]
    return out

## Function to get synonymous and non-synonymous diversity from a codon matrix
from utilities3.Data import genesitecount
gene_sites = genesitecount()
def snsdiv_wgscodonmat(gene,mat,gene_sites=gene_sites):
    # number of synonymous and non-synonymous sites
    S,N = gene_sites[gene][1:3]
    # no. of sequencing runs
    R = mat.shape[0]
    # no. of unique pairs
    P = binom(R,2)
    # get no. of synonymous and non-synonymous differences
    diffs = [ DnDs(mat[x,],mat[y,]) for x in range(0,R-1) for y in range(x+1,R)]
    # if any comparison fails OR no comparison was made, dont return diversity values
    if None in diffs or len(diffs) == 0: return
    Dn,Ds = [ sum(i) for i in zip(*diffs)]
    sdiv = round( Ds/(S*P), 4)
    nsdiv = round( Dn/(N*P), 4)
    return [gene,R,sdiv,nsdiv]

## Function to calculate symbol diversity ( general form of Nucleotide diversity)
# given a numpy matrix
def symbol_div(mat):
    # number of sites
    nsites = mat.shape[1]
    # indices of invariant columns
    rmcols = [ c for c in range(mat.shape[1]) if len(set(mat[:,c]))==1]
    # if all columns are invariant
    if len(rmcols) == mat.shape[1]:
        out = [0, 0]
    else:
        # remove invariant columns from the given matrix
        mat = np.delete( mat, rmcols,1)
        # no. of runs & no. of variant sites
        R,C = mat.shape
        # no. of pairwise comparison
        P = binom(R,2)
        S = sum( sum( mat[x,] != mat[y,]) for x in range(0,R-1) for y in range(x+1,R))
        D = round( S/(nsites*P), 4)
        out = [C,D] # C = variant sites, D = diversity
    return out

## Function to compute AA class diversity given a codon matrix
from utilities3.lists_dicts import mat_k2v
from utilities3.specific import mat_aa2class
from utilities3.Data import codon_table
ct = codon_table[11].forward_table
def aacdiv_cmat(gene,mat):
    # number of runs
    R = mat.shape[0]
    # convert codon matrix to amino acid matrix
    amat = mat_k2v(mat,ct)
    # get amino acid variant sites & diversity
    a_vs,a_div = symbol_div(amat)

    # convert amino acid matrix to AA class matrix
    acmat = mat_aa2class(amat)
    # get class variant sites & diversity
    ac_vs, ac_div = symbol_div(acmat)
    
    # output
    out = [ gene, R, a_vs, a_div, ac_vs, ac_div]
    return out

