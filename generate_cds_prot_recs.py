## generate multi-FASTA files for cds and protein sequence records of a gene across strains

# require
import os, sys, re
from utilities3.basic import *
from utilities3.Data import feature_table as ft
from Bio import SeqIO, SeqRecord
from joblib import Parallel, delayed

# start clock
start = time()

# arguments
nc = int(sys.argv[1]) # number of cores/threads

# paths
if not os.path.exists('cds'):
    os.mkdir('cds')
if not os.path.exists('prots'):
    os.mkdir('prots')

# load data
f_omat = 'ortho_matrix.tab'
omat = readcsv(f_omat,sep='\t')
# list of strains
strains = omat[0]
# and the matrix itself
omat = omat[1:]

# dict of reference protein ids with corresponding gene names
id_gname = { i[10]:i[14] for i in ft if i[0]=='CDS' and i[10] != ''}

## subset matrix to entries present in the reference
# reference assembly
ref = 'GCF_000005845.2_ASM584v2'
# column index of reference in the matrix
rx = strains.index(ref)
# indices of non-NA entries in the reference column
ixs = [ x for x,row in enumerate(omat) \
                if row[rx] != 'NA']
# get reference protein ids
ref_prots = [ omat[x][rx] for x in ixs]
# reduce matrix to above indices
omat = [ omat[x] for x in ixs]
# number of rows in the matrix
N = len(omat)

## Function to make a list of sequence records 
# note that this function doesn't make any checks; if the input files are not present, it throws an error
# else it continues and processes all records
def list_seqrecs(x,strain,matrix): 
                                    # x: column index in the ortholog matrix
                                    # strain: strain's genome's assembly accession
                                    # matrix: ortholog matrix
    print("%s: for %s\n"%(timer(start),strain))
    # extract corresponding column of protein ids
    ids = [ row[x] for row in matrix ]    
    # input FASTA file of coding sequences
    f_in = strain+'_cds_from_genomic.fna'
    # dict of sequence records
    seq_recs = SeqIO.to_dict( SeqIO.parse(f_in,'fasta'))
    # change above sequences ids to corresponding proteins ids
    seq_recs = { re.search('.*?_cds_?(.*)_[0-9]+',k).group(1):v for\
            k,v in seq_recs.items()}
    # records corresponding to ids
    out1 = [ seq_recs[i] if i in seq_recs.keys() else 'NA' for i in ids ] ## CDS

    # input FASTA file of protein sequences
    f_in = strain+'_protein.faa'
    # dict of sequence records
    seq_recs = SeqIO.index( f_in, 'fasta')
    # records corresponding to ids
    out2 =  [ seq_recs[i] if i in seq_recs.keys() else 'NA' for i in ids ] ## proteins
    return [out1,out2]

## Main
allrecords = Parallel( n_jobs = nc) ( delayed(list_seqrecs) (x,strain,omat)\
        for x,strain in enumerate(strains))
nucrecs = [ i[0] for i in allrecords]
protrecs = [ i[1] for i in allrecords]

# transpose above matrix, so that all records for a gene comes in the same list
tmat_n = [ [ j for j in i if isinstance(j,SeqRecord.SeqRecord) ] \
        for i in zip(*nucrecs)]
tmat_p = [ [ j for j in i if isinstance(j,SeqRecord.SeqRecord) ]\
        for i in zip(*protrecs)]

# for every gene
for x in range(N):
    # reference protein id
    prot = ref_prots[x]
    print("%s: for protein %i\n"%(timer(start),x+1))
    # corresponding gene name
    gene = id_gname[prot]
    # output file paths
    f_out_n = os.path.join( 'cds', gene+'.fasta')
    f_out_p = os.path.join( 'prots', gene+'.fasta')
    # write sequence records to a file
    SeqIO.write( tmat_n[x], f_out_n, 'fasta')
    SeqIO.write( tmat_p[x], f_out_p, 'fasta')

