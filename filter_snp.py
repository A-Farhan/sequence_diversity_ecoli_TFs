# extract and filter SNP data from relevant runs for selected genes

# require
import re, os, sys
from itertools import compress
from collections import Counter
from utilities3.basic import *
from utilities3.Data import genes_wo_additionalsites as genes
from utilities3.specific import *

# start clock
start = time()

# full path to the working directory
wdr = sys.argv[1]
os.chdir(wdr)

# file paths
fpa = 'genes_pa.tab'
fsnp1 = 'All_SNPs_annotated.txt'
fsnp2 = 'snps_filtered.tab'

# load data on gene presence/absence
pa_data = readcsv( fpa, sep = '\t')
# extract runs accessions from the header line
runs_pa = pa_data[0][1:]
# separate data from the header and convert binary sequence to int
pa_data = [ [i[0]]+[ int(j) for j in i[1:]] for i in pa_data[1:]]
# keep P/A data only for the above genes
pa_data = [ i for i in pa_data if i[0] in genes ]
print("%s: P/A data loaded from %s\n"%(timer(start),fpa))

# load SNP data
snp_data = readcsv( fsnp1, sep = '\t')
# extract header
snp_head = snp_data[0]
# separate SNP data from the header
snp_data = snp_data[1:]
# extract runs accessions from the header line
runs_snp = snp_head[2:-9]

# no. of runs
n_runs = len(runs_snp)
print("%s: This dataset has variants from %i runs\n"%(timer(start),n_runs))

# remove low quality variants
snp_data = rm_lq_snp( snp_data, n_runs)
print("%s: SNP data loaded and low quality variants removed, leaving %i snps\n"%(timer(start),len(snp_data)) )

# for each gene, identify the indices of runs 
# for which the gene was inferred present
# in the list of runs of SNP data
gene_runs = {}
# go through each row of gene P/A data
for n, row in enumerate(pa_data):
    # if not the first row and binary seq is same as the previous row
    if n > 0 and row[1:] == pa_data[n-1][1:]:     
        gene_runs[row[0]] = gene_runs[ pa_data[n-1][0]]
    else:
        gene_runs[row[0]] = get_runxs( inlist=runs_snp, maplist=runs_pa, binseq=row[1:])       
print( "%s: dict of runs indices generated for %i gene\n"%(timer(start),len(gene_runs)) ) 

# split SNP data by genes present in the genes P/A data
gene_snps = { k:[] for k in gene_runs.keys()}
# go through SNP data row by row
for row in snp_data:
    # if the gene is present in P/A data
    if row[-2] in gene_snps.keys():
                # and the variant is not upstream of the gene
        if row[-8] != 'upstream_gene_variant':
            # capture the row in the dict entry for the gene
            gene_snps[row[-2]].append(row)
print( "%s: SNP data was separated by %i genes\n"%(timer(start),len(gene_snps)) )

# extract and filter SNP data for genes in the P/A table
data_snp = []
# go through the above dict
for key, val in gene_snps.items():
    # no. of strains with the gene
    L = len( gene_runs[key])
    # go through each row of gene-specific data
    for row in val:
        # extract variant position
        pos = int( row[0].replace( 'NC_000913_', ''))

        # make sample sequence
        # use '.' for runs where gene is absent
        # Also include reference base
        S = ''.join( [ row[1]] + [ row[x+2] if x in gene_runs[key] \
                else '.' for x in range(n_runs) ])
        
        # find the most frequent base in the string
        B = Counter(S).most_common(1)[0][0] 
        # replace ambiguous bases with the most frequent base
        S = re.sub('[?]', B, S)
        
        # if the string has no variation other than ambiguous calls, skip
        set_chars = set(S)
        bases = list( filter( lambda x: x != '.', set_chars))
        if len(bases) <= 1: continue
        
        # remove reference base from the sample string
        S = S[1:]

        # add the following entry to the output table
        entry = [ key, pos, S, L] + row[-5:-3]
        data_snp.append(entry)

# order the above extracted data by position
data_snp.sort( key = lambda x: x[1])
print( "%s: SNP Data for genes in the P/A table was extracted and sorted by positions\n"%timer(start))
print("\t%i snp positions were extracted\n"%len(data_snp))

# write the extracted SNP data to file
writecsv( data = data_snp, fl = fsnp2, sep = '\t',\
        header = [ 'gene', 'pos', 'seq', 'presence', 'codon', 'aa'])
print( "%s: Extracted SNP data was written to %s\n"\
        %(timer(start),fsnp2))
