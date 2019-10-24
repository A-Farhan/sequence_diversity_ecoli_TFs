## generate protein MSAs and convert to codon MSAs
# require
import os, sys
from Bio.Align.Applications import ClustalOmegaCommandline \
        as clustalo
from joblib import Parallel, delayed
from subprocess import check_output
from Bio import SeqIO
from utilities3.basic import *

# start clock
start = time()

## arguments
try:
    nc = int( sys.argv[1])
except(IndexError):
    print("number of threads to use was not provided, default: 4")
    nc=4

# alignment program
aprog = sys.argv[2]
# program to convert amino acids alignments to codon alignments
cprog = sys.argv[3]

### Generate MSA ###
# paths
dr_n = 'cds'
dr_p = 'prots'
if not os.path.exists('pmsa'):
    os.mkdir('pmsa')

# list of cds FASTA files
fls_n = [ os.path.join( dr_n, i) for i in os.listdir(dr_n)]
# list of protein FASTA files
fls_p = [ os.path.join( dr_p, i) for i in os.listdir(dr_p)] 

## Function to run Clustal Omega parallely for multiple FASTA files
def clustal(infile,outfmt='fa',outdr='pmsa',prog=aprog):    
    # basename of input file path
    fname = os.path.basename(infile)
    # generate output file path
    outfile = os.path.join( outdr, fname)
    # if output already presents, then skip the computation
    if os.path.exists(outfile): return 
    # generate the command
    cmd = clustalo(prog,infile=infile,outfile=outfile,outfmt=outfmt)
    # run the command
    cmd() 
## run the clustalo function parallely
Parallel(n_jobs=nc)( map( delayed(clustal), fls_p) )
print('%s: MSA of proteins generated\n'%timer(start))

### Generate codon-alignments with protein MSA using PAL2NAL ###
dr1 = 'pmsa'
dr2 = 'cds'
dr3 = 'codon_msa'
if not os.path.exists(dr3):
    os.mkdir(dr3)

# list of files of protein MSA
fls = os.listdir(dr1)

# Function to convert protein MSA to CODON MSA
def prot2cds_msa(infile,dr1='pmsa',dr2='cds',dr3='codon_msa',prog=cprog):
    # output file path
    f_out = os.path.join( dr3, infile)
    # terminate, if output already exists
    if os.path.exists(f_out): return
    print("\t%s: For %s\n"%(timer(start),infile))

    # protein MSA file path
    fp = os.path.join( dr1, infile)
    precs = [ seq for seq in SeqIO.parse( fp, 'fasta')]
    # cds multi-FASTA file path
    fn = os.path.join( dr2, infile)
    nrecs = SeqIO.parse( fn, 'fasta')
    # ids of cds sequences
    ids = [ seq.id for seq in nrecs]
    # for every row in the MSA
    for y,rec in enumerate(precs):
        # generate id, name & description
        rec.id = ids[y]
        rec.name = ''
        rec.description = ''
    # create a temporary file to store the modified record
    SeqIO.write( precs, infile, 'fasta')
    # generate commands to convert amino acid alignments to codon alignments
    cmd = [ prog, infile, fn, '-output', 'fasta', '-codontable',\
            '11']
    # run command and capture output and write to file
    out = check_output(cmd).decode('ascii') 
    with open(f_out,'w') as flob:
        flob.write(out)
    # remove the temporary file
    os.remove(infile)
    
print("%s: Running PAL2NAL to generate codon alignments\n"%timer(start))
Parallel( n_jobs = nc) ( map( delayed(prot2cds_msa), fls)) 

