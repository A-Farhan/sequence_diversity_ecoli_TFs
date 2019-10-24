# Function to remove low quality variant calls from All SNP data
# requires header-removed data and # of sequencing runs
def rm_lq_snp( data, L):
    # index of column of basecode sequence
    ix = L+2
    cut = 0.1*L
    def main_fun( variant, ix=ix):
        S = variant[ix]
        return S.count('?') <= cut

    out = list( filter( main_fun, data))
    return out

# function to get list indices of sequencing runs present in SNP data
# if these are present in list of runs in which a gene is inferred "present"
from itertools import compress
def get_runxs( inlist, maplist, binseq):
    ixs = [ x for x, r in enumerate(inlist) if r in \
                list( compress( maplist, binseq))]
    return ixs

# rearrange a presence[1]/absence[0] matrix such that 
# 1. more filled rows appear first
# 2. identical rows appear together
def order_pamat( data):
    # cluster identical rows
    clusters = {}
    
    # go through each row from top to bottom
    for x in range(len(data)-1):
        
        # check if the index has already been captured in a cluster
        found = set( j for i in clusters.values() for j in i[1:])
        if x in found: continue
        
        # initiate a cluster(list) with the rowsum
        clusters[x]=[ sum(data[x][1:])]
        
        # make comparison with other rows excluding upstream rows
        for y in range(x,len(data)):
            if data[y][1:] == data[x][1:]:
                clusters[x].extend([y])

    # sort above cluster in decreasing order of their first value
    soclust = sorted( clusters.values(), reverse = True)

    # extract the order from sorted clusters
    order = [ j for i in soclust for j in i[1:]]

    # rearrage data in the order identified above
    out = [ data[i] for i in order ]
    return out

# test data
"""
data = [
        ['a',1,0],
        ['b',0,1],
        ['c',1,0],
        ['d',0,1]
    ]
"""

# recursive function to return a list of kmers
def get_kmers(k,A=['A','C','G','T']):
    if k == 1: return A
    out = [ i+j for i in get_kmers(k-1,A) for j in A]
    return out

# function to count all kmers in a given string
def count_kmers(dna, lkmer):
    kmer_counts = {} 
    kmers = [] 
    lseq = len(dna)
    # generate all kmers present in the sequence
    for start in range(lseq-lkmer+1):
        kmer = dna[start:start+lkmer]
        kmers.append(kmer)
    for a_kmer in kmers:
        old_count = kmer_counts.get(a_kmer,0)
        new_count = old_count + 1
        kmer_counts[a_kmer] = new_count
    return kmer_counts
 
# function to generate N shuffles of a sequence preserving 
# k-mer frequencies
def shuffle_w_kmerf(p,k,N=1):
    # length of promoter
    l = len(p)
    # get kmer frequencies in the given seq
    kf = count_kmers( p, k)
    # convert above to a list of kmers repeated a/c to above freqs
    ks = [ k for k, v in kf.items() for i in range(v)]
    # initialize list to hold shuffled sequences
    shufseqs = []
    # repeat until N sequences are generated
    while len(shufseqs) < N:        
        # randomly pick 1 kmer from the above list
        k1 = random.sample( ks, 1)[0]
        # initialize a sequence with it
        seq = k1
        # repeat until the seq is as long as the original promoter
        while len(seq) < l:
            # limit selection of next kmer to those starting with 
            # last character of the seq
            lim_ks = [ i for i in ks if i[0]==seq[-1]]
            kn = random.sample( lim_ks, 1)[0]
            # append the kmer to sequence excluding the first char
            seq = seq + kn[1:]
        
        # capture a unique permutation
        if seq not in shufseqs:
            shufseqs.append(seq)
    return shufseqs

## function to convert a list of nucleotides to a list of codons
n2c = lambda seq: [ ''.join(seq[x:x+3]) for x in range(0,len(seq),3)]

## function to convert a list of codons to a list of nucleotides
c2n = lambda seq: [ j for i in seq for j in i]

### Function to get site counts for a given sequence record
from Bio import AlignIO as align
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table
from utilities3.Data import codonsitecount
site_dict = codonsitecount()
def recsitecount(seqrec,ct=codon_table[11],nucs='ACGT',site_dict=site_dict):
    # check for ambiguous letters
    check = any( i not in nucs for i in seqrec.seq if i != '-')      
    if check: return
    
    # nucleotide sequence
    nseq = str(seqrec.seq)
    # length of sequence
    L = len(nseq)
    # codon sequence
    cseq = n2c(nseq)
        
    # check for presence of internal stop codons
    check = any( i in ct.stop_codons for i in cseq)
    if check: return

    # genome id
    Id = seqrec.id.split('_cds_')[0][4:]
        
    # sequence without gaps
    cseq = [ i for i in cseq if '-' not in i]
    # sites without gaps
    T = len(cseq)*3
    # sites for all codons
    sites = [ site_dict[codon] for codon in cseq]
    S,N,F = ( round( sum(i), 2) for i in zip(*sites))
    out = [Id,L,T,S,N,F]
    return out

## Function to generate a dict of sequencing runs with their genes
# given a P/A matrix
import os, re
from itertools import compress
from utilities3.basic import * 
def rungenes_dict(dr):
    f_gpa = os.path.join( dr, 'genes_pa.tab')
    # matrix with runs on rows
    mat = [ list(i) for i in zip(*readcsv(f_gpa,sep='\t'))]
    # list of genes
    genes = mat[0][1:]
    # no. of runs
    N = len(mat)-1
    out = { mat[x+1][0]: list(compress(genes, map( int, mat[x+1][1:])))\
        for x in range(N)}
    return out

## Function(recursive) to get count of Non-synonymous and Synonymous 
# differences for a single codon
from itertools import permutations
from math import factorial
def cDnDs(c1,c2,ct=codon_table[11]):
    table = ct.forward_table
    stops = ct.stop_codons
    """
    # check its a codon
    assert len(c1) == len(c2)
    assert len(c1) == 3
    """
    # initialize output variables
    Dn,Ds = 0,0

    # check if any given codon is a stop codon
    if c1 in stops or c2 in stops: return None
    
    # get differing positions
    dp = [ y for y in range(3) if c1[y] != c2[y]] 
    # count no. of differences
    nd = len(dp)

    # if no differences then both Dn & Ds output is zero
    if nd == 0: 
        return [Dn,Ds]
    
    # if only one difference, then it's either Dn or Ds is 1
    if nd == 1:
        # corresponding amino acids
        a1 = table[c1]
        a2 = table[c2]
        # Non-synonymous if amino acids are not the same
        if a1 != a2: Dn = 1
        else: Ds = 1
        return [Dn,Ds]
    
    ## For more than 1 difference
    
    # no. of permutations
    np = factorial(nd)
    # list of permutations of above positions
    perms = permutations(dp)
    # initialize list to hold mutations in paths
    pmut = [ [] for x in range(np)]
    # for every permutation of mutated positions
    for t,tup in enumerate(perms):
        # initialize path with starting codon
        path = [c1]
        # initialize value of intermediate codon 
        # with starting codon
        ic = c1
        # for every position, generate an intermediate codon
        for p in tup:
            ic = ''.join( [ c2[x] if x == p else ic[x] \
                    for x in range(3)])
            # and add to path
            path.append(ic)
        # go through path, 2 codons at a time, get their difference
        for x in range(nd):
            res = cDnDs(path[x],path[x+1]) # path length is 1 greater than # mutations
            pmut[t].append(res)  
            # if at any stage, a stop codon is being compared
                # don't proceed on this path, move on to next
            if res is None: break
            # normally
            Dn += res[0]
            Ds += res[1]
    # exclude paths with stop codons
    pmut = [ i for i in pmut if None not in i]
    # turn above into an unnested list
    allmut = [ j for i in pmut for j in i]
    # get average Dn,Ds values 
    res = [ round( sum(i)/len(pmut), 2) for i in zip(*allmut)] 
    return res

# Function to get Non-synonymous and Synonymous differences
# given two coding sequences
def DnDs(cs1,cs2):
   # no. of codons
    nc = len(cs1)
    # initialize variable to hold count of differing sites
    Dn,Ds = 0,0
    # move through sequences 1 codon at a time
    for x in range(nc):
        # codons
        c1 = cs1[x]
        c2 = cs2[x]
        counts = cDnDs(c1,c2)
        Dn += counts[0]
        Ds += counts[1]
    return [Dn,Ds]

###########################
# Find minimum distance threshold to decluster, given a distance matrix of objects
# require
import numpy as np
import pandas as pd
from utilities3 import basic
def decluster_distance_threshold(distmat,minslope):
    # rename arguments
    dm = distmat
    b = minslope
    # no. of objects in the distance matrix
    N = len(dm)
    # sorted list of pairwise distances
    dlist = sorted( dm[r,c] for r in range(N-1) for c in range(r+1,N))
    # no. of non-redundant pairwise comparison
    L = len(dlist)
    # log CDF of above distance values
    lcf = [ basic.log( sum( map( lambda x: x <= i, dlist))/L) \
            for i in dlist]
    lcf = list( filter( None, lcf))
    # convert to dataframe with distance & log CDF as columns
    df = pd.concat([ pd.DataFrame(dlist), pd.DataFrame(lcf)], axis=1) 
    df.columns = ['dist','LCF']
    # keep rows corresponding to first instance of each distance
    df = df[ ~df.duplicated('dist')]
    # covert to numpy array
    mat = df.values
    # remove rows corresponding to zero distance
    mat = mat[ mat[:,0] != 0, ]
    # no. of rows in above matrix
    R = len(mat)
    # initialize list of slopes
    slopes = [ 0 for _ in range(R)]
    slopes[0] = mat[0,1]/mat[0,0]
    slopes[1:] = [ (mat[r,1]-mat[r-1,1]) / (mat[r,0]-mat[r-1,0]) \
            for r in range(2,R)]
    # get indices of slopes that satisfy given threshold
    ixs = [ x for x,i in enumerate(slopes) if i < b and i >= 0]
    # get corresponding distance
    out = mat[ ixs[0], 0]
    return out


## Function to convert amino acid matrix to aa class matrix
import numpy as np
from utilities3.Data import aa_class as aac
aac = { k[2]:v for k,v in aac.items()}
def mat_aa2class(amat,aac=aac):
    a = amat
    a_flat = np.ndarray.flatten(a)
    b = [ [k for k,v in aac.items() if x in v][0] if x != '*' else\
            'x' for x in a_flat]
    c = np.reshape( b, a.shape)
    return c

