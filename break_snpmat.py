# generate separate mutations tables for each run from snp matrix
# the tables will have columns: gene, codon_num(1-based), alt

## variants table are not generated for sequencing runs without gene lists

# require
import os, sys, re, operator
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from joblib import Parallel,delayed
from utilities3.basic import * 
from utilities3.wgs_diversity import rungenedict 
from utilities3.Data import feature_table as ft, \
        genes_wo_additionalsites as genes

# start clock
start = time()

# arguments
# number of threads
nt = int(sys.argv[1])
# full path to the working directory
wdr = sys.argv[2]

# paths
fsnp1 = os.path.join( wdr, 'All_SNPs_annotated.txt')
fsnp2 = os.path.join( wdr, 'snps_filtered.tab')

# list of sequencing runs, extracted from the SNP matrix header
seqruns = [ i for i in readcsv( fsnp1, sep='\t')[0] if 'RR' in i]
# curated SNP matrix data
snpmat = readcsv( fsnp2, sep='\t')[1:]

## code to extract codon position
# compile regular expression
rc = re.compile('p.[A-Z][a-z]{2}([0-9]+)')
def extract_codonpos(String):
    M = rc.match(String)
    if M is None: return
    else: return M.group(1)

## code to fix orientation of sample sequences
# complementing for genes on the negative strand
gene_strand = { i[14]:i[9] for i in ft if i[14] in genes}
def fix_orientation(sampleseq,gene,gene_strand=gene_strand):
    if gene_strand[gene] == '+':
        return sampleseq
    else:
        return str(Seq( sampleseq,IUPAC.unambiguous_dna).complement())

## code to process snp matrix
# keep gene name, codon position, reference codon, variants' sequence
def process_snpmat(row,genes=genes): # row refers to a row of filtered SNP matrix
    # skip the row if gene not in the selected set
    if row[0] not in genes: return

    # codon position
    codonpos = extract_codonpos(row[-1]) # the last column in the SNP matrix has the annotation 
    # fix orientation of the sample sequence
    sampleseq = fix_orientation(row[2],row[0])
    # reference codon
    ref = row[-2].split('/')[0] # the first part of the second last column in the SNP matrix
                                # the variant nucleotide is capitalized
    return [ ':'.join( [row[0],codonpos]),ref,sampleseq] # join gene & codon position 
# apply above function to all rows of curated SNP matrix parallely
processed = filter( None, Parallel(n_jobs=nt) \
        ( map( delayed(process_snpmat), snpmat)))    
print("%s: SNP matrix was processed\n"%timer(start))

# split the data by first row ( gene + codon position)
pos_snp = split_data(processed,0)
print("%s: the SNP data was split by gene & codon position\n"%timer(start))

## code to convert nucleotide sample sequence to codon sequence
# joining variant rows corresponding to the same codon
# Function
def nseq2cseq(nucrows): # nucrows refers to rows of SNP matrix split by gene and codon position
                        # st. rows falling in the same codon are supplied together as nested lists
    for x,i in enumerate(nucrows):
        # from the first row, extract gene,position and reference codon ( can be converted to all capital letters)
        if x==0: 
            gene,pos = i[0].split(':')
            ref = i[1].upper()

        # mutated position in the reference codon, identified by the sole capitalized letter
        mp = [ y for y,j in enumerate(i[1]) if j.isupper()][0]
        
        # generate variant codon sample sequences
        # sequentially, starting from the reference
        if x==0:
            alts = [ ''.join([ n if y == mp else ref[y] for \
                    y in range(3)]) for n in i[2]]
        else:
            alts = [ ''.join([ n if y == mp else alts[z][y] for \
                    y in range(3)]) for z,n in enumerate(i[2])]
    # replace all codons with gaps by '.'
    alts = [ '.' if '.' in j else j for j in alts]
    return [gene,int(pos),ref] + alts # converting position to integer is required, for a sorting done later

# Apply above function to dict of rows corresponding to the same codon parallely
codonseqs = Parallel(n_jobs = nt) \
        ( map( delayed(nseq2cseq), pos_snp.values()))

# sort above by the gene names and codon position
codonseqs.sort(key = operator.itemgetter(0,1))
print("%s: codon sequences were generated\n"%timer(start))

## code to separate mutations for each sequencing run
run_genes = rungenedict(seqruns) # dict of sequencing runs and their genes
# Function
def seqrun_vars(index,name,codonseqs=codonseqs,run_genes=run_genes): 
                        # index and name refers to a sequencing run in the list of runs
                        # codonseqs refers to the matrix of sample sequences in the codon form (generated above)
    # output file path for the sequencing run
    f_out = os.path.join( name, 'variants.tab')
    # Exit, if file exists
    if os.path.exists(f_out): return
    # Exit if the run absent from dict of runs with genes
    if name not in run_genes.keys(): return 
    
    # run's index in the list of runs
    x = index
    # runs' accession(name)
    sr = name
    # initialize output list
    out = []
    # for every row in the codonseq matrix
    for row in codonseqs:
        # reference codon
        ref = row[2]
        # codon in the sequencing run
        alt = row[x+3]
        # skip the row if run's codon is same as reference
        if alt == ref: continue
        # skip the row if gene is not present in the run
        if row[0] not in run_genes[sr]: continue
        # otherwise, extract features from the row
        out.append(row[:2]+[alt])
    writecsv(data=out,fl=f_out,sep='\t',\
            header=['gene','codonpos','var'])
# Apply above function to the list of sequencing runs parallely
Parallel(n_jobs=nt) ( delayed(seqrun_vars) (x,i) \
                    for x,i in enumerate(seqruns))
print("%s: separate files generated for all sequencing runs\n"%timer(start))
