## Selection of refseq strains to avoid redundancy ##

# require
import os, sys
from subprocess import call
from utilities3.basic import *
from utilities3.lists_dicts import *

# start clock
start = time()

# arguments
db_dr = sys.argv[1] # directory where genomes are stored

# paths
f_st = 'eck12mg1655_machtman_mlst'

### Extract sequences of 7 genes corresponding to E coli K-12 MG1655 strain's ST type based on Mark Achtman's MLST scheme
# list of gene names and allele number for reference e coli in this scheme
alnum = readcsv(f_st)
# No. of mlst genes
N = len(alnum)
# their names
genes = [ i[0] for i in alnum]
# do the following if reference allele sequence file doesn't exist
f_ref = 'ref_mlst_gfrag.fas'
if not os.path.exists(f_ref):
	# dict of genes with sequences
	gene_seqs = {}
        # for every gene
	for item in alnum:
	    f_in = item[0]+'.fas'
	    gene_seqs[item[0]] = multifasta_dict(f_in)[''.join(item).upper()]
	# write to file
	with open(f_ref,'w') as flob:
	    for gene,seq in gene_seqs.items():
	        entry = '>'+gene+'\n'+seq+'\n'
	        flob.write(entry)

## length of reference alleles
gene_length = { k:len(v) for k,v in multifasta_dict(f_ref).items()}

### Get percent identity of MLST fragments of refseq strains to the reference

## perform BLASTn
# cut-offs
ev = 1e-05 # evalue
pid = 70 # percent identity
po = 90 # percent overlap

# list of refseq strains
fstrain = 'refseq_strains.txt'
with open(fstrain) as flob:
    strains = [ i.strip('\n') for i in flob]

# directory to store blastn output
dr_out = 'blastn_mlst'
if not os.path.exists(dr_out):
    os.mkdir(dr_out)

for x,s in enumerate(strains):
    sub_f = os.path.join( db_dr, s + '_genomic.fna')
    f_out = os.path.join( dr_out, s + '.tab')
    if not os.path.exists(f_out):
        print("%s: Running BLASTn against %s\n"%(timer(start),sub_f))
        cmd = ['blastn', '-query', f_ref, '-subject', sub_f, '-outfmt', '6', '-out', f_out, \
                '-max_target_seqs', '1', '-max_hsps', '1', '-evalue', str(ev), \
                '-perc_identity', str(pid), '-qcov_hsp_perc', str(po)]
        call(cmd)

## generate a matrix of percent identities
# file to store matrix
f_mat = 'mlst_identity_matrix.tab'
# list of blastn output files
fls = os.listdir(dr_out)
# no. of assemblies
A = len(fls)
if not os.path.exists(f_mat):
	# initialize the matrix
	idenmat = [ [ 'NA' for c in range(N+1)] for r in range(A)]
	# fill the matrix
	for r,f in enumerate(fls):
	    asm = f.replace('.tab','')
	    f_in = os.path.join( dr_out, f)
	    gene_iden = { i[0]:[ i[2], \
	            str( float(i[3])/gene_length[i[0]]*100)]\
	            for i in readcsv(f_in,sep='\t')}
	    idenmat[r][0] = asm
	    for c,gene in enumerate(genes):
	        if gene in gene_iden.keys():
	            entry = ','.join( [ i.split('.')[0] \
	                    for i in gene_iden[gene]])
	            idenmat[r][c+1] = entry
	# write above matrix to a file
	writecsv(data=idenmat,fl=f_mat,sep='\t',\
	        header=['assembly']+genes)


### select strains based on the above matrix such that
# no two strains have a difference from the reference
# of less than 1% in percent identity and overlap of all MLST genes

# input file
f_in = f_mat
# output file
f_out = 'selected_refseq_strains.txt'
# reference assembly accession
ref = 'GCF_000005845.2_ASM584v2'

# load identity matrix
idenmat = readcsv(f_in,sep='\t')[1:]

# Do the following, ONLY if the strain set file is not already present
if not os.path.exists(f_out):
	# remove the strains with NA
	idenmat = [ i for i in idenmat if 'NA' not in i]
	# convert above list to a dict
	asm_data = { i[0]:i[1:] for i in idenmat}
	# find unique seqeunces of identities
	unq = dedup( asm_data.values())
	# club assemblies by above sequences
	asmsets = []
	for i in unq:
	    entry = []
	    for k,v in asm_data.items():
	        if v==i: entry.append(k)
	    asmsets.append(entry)
	
	# select one strain from each set
	strains = [ i[0] if ref not in i else ref for i in asmsets]
	with open(f_out,'w') as flob:
	    for i in strains: flob.write(i+'\n')
