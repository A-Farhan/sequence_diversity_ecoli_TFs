### Argument parsing ###
import argparse
parser = argparse.ArgumentParser(description="Generates a Matrix of genes \
        with orthologs' ids present in the given set of genomes. Singletons are not included. \
        Before running the program, make sure that all <genome>_protein.faa and <genome>_feature_table.txt files\
        are present in the current working directory,either copied or linked. \
        The script should be called from within the working directory.")
parser.add_argument('-P', '--program', help='full path to the directory with BLASTp program,\
                    including the trailing "/"', default='')
parser.add_argument('-C', '--num_threads', help='number of threads to use',default='1')
mandate_arg = parser.add_argument_group('mandatory arguments')
mandate_arg.add_argument('-i', '--input_file', help="file with a list of strains' assembly accessions",required=True)
blast_arg = parser.add_argument_group('BLASTp parameters')
blast_arg.add_argument('--pid', help='minimum percent identity (default: %(default)s)',default='50')
blast_arg.add_argument('--po', help='minimum query coverage (default: %(default)s)',default='50')
blast_arg.add_argument('--ev', help='maximum e-value (default: %(default)s)',default='1e-05')
blast_arg.add_argument('--matrix', help='substitution scoring matrix (default: %(default)s)',default='BLOSUM80')

args = parser.parse_args()
######################################

# require
import os, sys
from subprocess import call
from collections import Counter
from joblib import Parallel, delayed
from utilities3.basic import *
from utilities3.lists_dicts import *
from orthogene_groups3.functions import *

# start clock
start = time()

# arguments
f_strains = args.input_file
path_bl = args.program
nc = args.num_threads
spacer='\t'

# blastp parameters
ev = args.ev
pid = float(args.pid)
po = float(args.po)
matrix = args.matrix

# list of strains
with open(f_strains) as flob:
    strains = [ i.strip('\n') for i in flob]

# no. of strains
N = len(strains)
print("%s: Input data has %i strains\n"%(timer(start),N))

# directory to store local databases built from protein sequences 
ldb = 'local_db'
if not os.path.exists(ldb):
    os.mkdir(ldb)
# directory to store reciprocal genome BLASTp output
bdr = 'pwgblast'
if not os.path.exists(bdr):
    os.mkdir(bdr)
# directory with pairwise orthologs
dr_pwo = 'pw_ortho'
if not os.path.exists(dr_pwo):
    os.mkdir(dr_pwo)
# directory for ortholog matrices
dr_om = 'ortho_matrices'
if not os.path.exists(dr_om):
    os.mkdir(dr_om)

###################################################################

## Make BLASTp local databases ##
Parallel(n_jobs = int(nc)) ( delayed(make_blastdb) (strain=i,dr=ldb,path_bl=path_bl) \
        for i in  strains)

## Perform pairwise reciprocal BLASTp ##
for x in range(N):
    n1 = strains[x]
    # blastp query file
    f_query = n1+'_protein.faa'
    for y in range(N):
        n2 = strains[y]
        if n1 == n2: continue
        # blastp database name
        dbp = os.path.join( ldb, n2)
        outf = os.path.join( bdr, n1+'__'+n2+'.tab')
        
        # move to next, if BLASTp output already exists
        if os.path.exists(outf): continue
        
        # ~ 26 sec per comparison
        print("%s: BLASTp %s with %s\n"%(timer(start),n1,n2))
        cmd = [path_bl+'blastp', '-query', f_query, '-db', dbp,\
            '-outfmt', '6', '-out', outf, '-max_target_seqs', '1',\
            '-max_hsps', '1', '-evalue', ev, \
            '-num_threads', nc, '-matrix', matrix]
        call(cmd)

## Find Bidrectional Best Hits ##

# Function constructor for filtering by overlap
def makefilter(list1,list2,po=po):
    def po_filter(item):
        k1,k2 = item[:2]
        # if k1 not in list1 or k2 not in list2: 
        #     return False
        overlap = float(item[3])   
        maxl = max( list1[k1], list2[k2])
        percent = overlap/maxl * 100
        if percent >= po:
            return True
        else:
            return False
    return po_filter

# for every pair of strains
for x in range(N-1):
    # strain 1
    n1 = strains[x]
    # feature table 
    f_ft1 = n1+'_feature_table.txt'
    ft1 = filter( lambda x: x[0]=='CDS' and len(x[-2]) > 0, \
            readcsv(f_ft1,sep='\t'))
    # dict of protein's lengths
    pl1 = { i[10]:int(i[-2]) for i in ft1}
    
    for y in range(x+1,N):
        # strain 2
        n2 = strains[y]
        # move to next, if output present
        outf = os.path.join( dr_pwo, n1+'__'+n2+'.tab')
        if os.path.exists(outf): continue
        
        print("%s: Extracting orthologs between %s and %s\n"%(\
                timer(start),n1,n2))
        # feature table 
        f_ft2 = n2+'_feature_table.txt'
        ft2 = filter( lambda x: x[0]=='CDS' and len(x[-2]) > 0,\
                readcsv(f_ft2,sep='\t'))
        # dict of protein's lengths
        pl2 = { i[10]:int(i[-2]) for i in ft2}

        # forward and reverse BLASTp result paths
        f_fw = os.path.join( bdr, n1+'__'+n2+'.tab')
        f_rev = os.path.join( bdr, n2+'__'+n1+'.tab')
        # corresponding data
        d_fw = readcsv(f_fw,sep='\t')
        d_rev = readcsv(f_rev,sep='\t')
        
        # filter based on identity
        fw1 = [ i for i in d_fw if float(i[2]) >= pid]
        rev1 = [ i for i in d_rev if float(i[2]) >= pid]
        # filter based on overlap
        filter_fw = makefilter( pl1, pl2)
        filter_rev = makefilter( pl2, pl1)
        fw2 = filter( filter_fw, fw1)
        rev2 = filter( filter_rev,rev1)
        # list of first two columns joined
        fw = [ i[0]+'__'+i[1] for i in fw2]
        rev = [ i[1]+'__'+i[0] for i in rev2]

        # keep only those in forward which are present in reverse
        bd = [ i.split('__') for i in fw if i in rev]
        # write to file
        writecsv(data=bd,fl=outf,sep='\t')

## Make query-based orthologues matrices ##
print("%s: Generating query-based orthologous matrices\n"\
        %timer(start))
Parallel(n_jobs = int(nc)) \
    ( delayed(query_orthomat) (query_x=x,inlist=strains,outdr=dr_om,dr_pwo=dr_pwo) for x in range(N) ) 

print("%s: Loading orthologous matrices\n"%timer(start))
holder = []
for x in range(N):
    f_in = os.path.join(dr_om, strains[x]+'.tab')
    mat = readcsv(f_in,sep='\t')
    # exclude proteins only present in the query
    holder.append( [ row for row in mat if not \
            row.count('NA') == N-1])

## Find the strict core ##
f_strictcore = 'strict_core.tab'
if not os.path.exists(f_strictcore):
    print("%s: Generating  strict core\n"%timer(start))
    strict_core = [ i for i in holder[0] if all( [ i in holder[x] \
                for x in range(1,N)])]
    writecsv(fl=f_strictcore,data=strict_core,\
            sep='\t',header=strains)
else:
    print("%s: strict core file already present. Loading\n"%timer(start))
    strict_core = readcsv(f_strictcore,sep='\t')[1:]
print("%s: Excluding strict core rows from all matrices\n"%timer(start))
holder = [ [ j for j in i if j not in strict_core] for i in holder ]

## Find closed groups ##

# every member find all other members
f_close = 'closed_groups.tab'
if not os.path.exists(f_close):
    print("%s:Initialising list to hold closed groups\n"%timer(start))
    closeg = []
    rejects = []
    # for every query matrix
    for x in range(N-1):
        # for every row of this matrix
        print("%s: Looking for closed groups in genome %i\n"%(spacer,x+1))
        for row in holder[x]:
            # skip the row if already included or rejected
            if row in closeg or row in rejects: continue
            
            # indices of the matrices where row reports a hit
            ixs = [ y for y in range(N) if y != x and row[y] != 'NA']
            # check if row is NOT found in above indices
            not_close = any( row not in holder[y] for y in ixs)
            # if not, reject the row
            if not_close:
                rejects.append(row)
            # otherwise, include
            else:
                closeg.append(row)
    # write to existing file   
    writecsv(data=closeg,fl=f_close,sep='\t',header=strains)
else:
    print("%s: Closed group file already exists. Loading\n"%timer(start))
    closeg = readcsv(f_close,sep='\t')[1:] 

print("%s: Excluding closed groups\n"%timer(start))
holder = [ [ j for j in i if j not in closeg] for i in holder]

## Extract orthologs from open groups ##
# those where not all members find all other members, but only the closed part is retained
f_open = 'open_groups.tab'
if not os.path.exists(f_open):
    openg = find_opengroups(holder)
    # write matrix to file
    writecsv(data=openg,fl=f_open,sep='\t',header=strains)
else:
    print("%s:%s already exists\n"%(timer(start),f_open))
    openg = readcsv(f_open,sep='\t')[1:]

# generate complete orthologous matrix
f_out = 'ortho_matrix.tab'
if not os.path.exists(f_out):
    mixed = readcsv(f_open,sep='\t')[1:]
    out = strict_core + closeg + openg
    writecsv( data=out, fl=f_out, sep = '\t', header = strains)
print("%s: complete orthologous matrix generated\n"%timer(start))
