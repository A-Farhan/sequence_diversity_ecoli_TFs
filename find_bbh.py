# find BBH of a given protein across a list of proteomes
# and create a multi-FASTA file

## Note: Run this script from within a directory dedicated to Across species comparison 
## make sure the input gene's protein sequence from the reference genome is present in the directory in FASTA format

# require
import sys, re
from subprocess import call
from utilities3.basic import *
from Bio import SeqIO
# arguments
gene = sys.argv[1] # input gene name
ev = 1e-06

# start clock
start = time()

## paths
drp = 'refprogam'
dint = 'intfiles' # directory to dump intermediate files
fgs = gene + '_prot.fasta'
fo_tab = gene + '_hits.tab'
fo_mf = gene + '_bbh.fasta'

# reference proteome 
ref = 'UP000000625_83333'
# file path
fref = os.path.join( drp, ref+'.fasta')
# protein id in the reference proteome
id_ref = [ k for k,v in SeqIO.index(fref,'fasta').items() \
            if bool(re.search('GN='+gene,v.description))][0]

# list of proteomes present, excluding reference proteome
pomes = [ i.replace('.fasta','') for i in os.listdir(drp)]
pomes.remove(ref)

# initialize list of hits
hits = [ [ref,id_ref] ]
## run PHMMER and find BBH
for c,pome in enumerate(pomes):
    print("%s: For proteome %i\n"%(timer(start),c+1))
    # input file for proteome
    fin = os.path.join( drp, pome+'.fasta')
    # output file for phmmer forward
    fo_pf = os.path.join( dint, '__'.join([gene,pome]) )
    # run phmmer, if this file doesn't exist
    if not os.path.exists(fo_pf):
        cmd = ['phmmer','-E',str(ev),'--tblout',fo_pf,fgs,fin]
        call(cmd)
    # load results for phmmer forward, if no results then skip
    res = readcsv(fo_pf,sep=' ')
    if len(res) == 0: continue
    # best hit
    bh = res[0][0]
    # truncated name
    bh_trunc = bh.split('|')[-1]
    # output fasta file for the best hit
    ffbh = os.path.join( dint, bh_trunc + '.fasta')
    # write the bh's sequence to a separate fasta file,
    # if it doesn't already exist
    if not os.path.exists(ffbh):
        # sequence record for the above hit
        bh_seqrec = SeqIO.index(fin,format='fasta')[bh]
        with open(ffbh,'w') as flob:
            SeqIO.write(bh_seqrec,flob,'fasta')
    # output file for phmmer reverse
    fo_pr = os.path.join( dint, '__'.join([pome,gene]) )
    # run phmmer, if this file doesn't exist
    if not os.path.exists(fo_pr):
        cmd = ['phmmer','-E',str(ev),'--tblout',fo_pr,ffbh,fref] 
        call(cmd)
    # load results for phmmer forward, if no results then skip
    res = readcsv(fo_pr,sep=' ')
    if len(res) == 0: continue
    # reverse best hit
    rbh = res[0][0]
    # if the above hit is same as the original query
    if rbh == id_ref:
        hits.append([pome,bh])

# write table of hits to file
writecsv(data=hits,fl=fo_tab,sep='\t',\
        header=['proteome','protein'])

## create a multi-FASTA file of hits
hits_seq = []
for x,row in enumerate(hits):
    pome = row[0]
    prot = row[1]
    fin = os.path.join( drp, pome+'.fasta')
    hitrec = SeqIO.index(fin,format='fasta')[prot]
    hits_seq.append(hitrec)
with open(fo_mf,'w') as flob:
    SeqIO.write(hits_seq,flob,'fasta')

