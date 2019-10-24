import os, re
from utilities3 import basic
path = os.path.dirname(basic.__file__)
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table

## reference feature table
f_in = os.path.join( path, 'reference.ft')
data = basic.readcsv(f_in,sep='\t')[1:]
feature_table = [ i for i in data if i[0]=='CDS' and i[-1]!='pseudo']

## list of reference genes with no additional site to the reference genome
# the sequences include the stop codon
f_in = os.path.join( path, 'noadditionalsitegenes.txt')
with open(f_in) as flob:
    genes_wo_additionalsites = [ i.strip('\n') for i in flob]

## list of genes with same site in all refseq strain as the reference 
# the sequences include the stop codon
f_in = os.path.join( path, 'samesitegenes.txt')
with open(f_in) as flob:
    genes_w_samesites = [ i.strip('\n') for i in flob]

## list of core genes from refseq strains with high nucleotide diversity
f_in = os.path.join( path, 'highdivgenes.txt')
with open(f_in) as flob:
    refhdgenes = [ i.strip('\n') for i in flob]

## list of core genes from wgs strains with high nucleotide diversity
f_in = os.path.join( path, 'raw_highdivgenes.txt')
with open(f_in) as flob:
    rawhdgenes = [ i.strip('\n') for i in flob]

## dict of sequence records of genes from the reference genome
f_in = os.path.join( path, 'reference.fna')
seqrecs = SeqIO.parse(f_in,'fasta',alphabet=IUPAC.unambiguous_dna)
# regular expression object to extract gene name
rc = re.compile('\[gene=([a-z]{3}[A-Z]?)\]')
# function using above object to extract gene name from the header
def extractgenename(header):
    m = rc.search(header)
    if m is None: return 
    else: gname = m.group(1) 
    return gname
seqrecord = {}
for rec in seqrecs:
    gname = extractgenename(rec.description)
    if gname is not None: 
        seqrecord[gname] = rec.seq

## Synonymous and Non-synonymous site count for all codons according to CodonTable(Bacteria)
def codonsitecount(ct = codon_table[11],nucs = 'ACGT',ixs=range(3)):
    table = ct.forward_table
    stops = ct.stop_codons
    codons = table.keys()
    out = { k:() for k in codons}
    for codon in codons:
        # corresponding amino acid
        aa = table[codon]

        # Get Synonymous Site Count
        # initialize list to hold fraction of synonymous changes at each site of a codon with zeros
        f = [ 0 for _ in ixs]
        # for every site in the codon
        for i in ixs:
            # produce all possible alternate codons 1 mutation away
            altc = [ ''.join( [ n if x == i else codon[x] \
                        for x in ixs]) for n in nucs if n != codon[i]]
            # corresponding amino acids for above codons
            alta = [ table[c] for c in altc if c not in stops]
            # fraction of above codons coding for same amino acid
            f[i] = float(alta.count(aa))/len(alta)
        S = round(sum(f),2)

        # Get non-synonymous site count
        N = 3 - S
        
        # Get count of 4-fold degenerate sites
        # alternate codons with 1 mutation at the last site
        altc = [ codon[:2]+n for n in nucs \
                        if n != codon[2]]
        # corresponding amino acids for above codons
        alta = [ table[c] for c in altc if c not in stops]
        # check if all 3 codons code for the same amino acid as original
        if alta.count(aa) == 3: 
            F = 1
        else:
            F = 0
        out[codon] = (S,N,F)
    return out

## sites count for all reference genes
def genesitecount(rec=seqrecord,sites=codonsitecount()):
    out = {}
    for k,v in rec.items():
        # nucleotide sequence
        nseq = str(v)[:-3]
        # length of the sequence OR total no. of sites
        L = len(nseq)
        # corresponding codon sequence
        cseq = [ nseq[x:x+3] for x in range(0,L,3)]
        # list of tuples of site counts for each codon
        counts = [ sites[codon] for codon in cseq]
        # transpose and sum to get total site count of each type
        (S,N,F) = ( round(sum(i),2) for i in zip(*counts))
        # recalculate N by subtracting S from total sites
        # to fix decimal errors
        N = round(L-S,2)
        out[k] = [L,S,N,F]
    return out

### Code to make lists of TFs and target genes ###
def reg_genes():
    # paths
    fnet = os.path.join( path, 'network_tf_tu.txt')
    ftu = os.path.join( path, 'TUSet.txt')
    ftf = os.path.join( path, 'tf_exp_regulon.txt')

    # load data
    net = basic.readcsv(fnet,sep='\t')
    net = [ [ i[0], i[1].split('[')[0]]  for i in net]

    tus = basic.readcsv(ftu,sep='\t')
    # convert above to dict of TU name and constituent genes
    tus = { i[1]:i[3].split(',') for i in tus}

    exptf = basic.readcsv(ftf,sep='\t')
    # convert above to dict of TF name with gene names
    tf_gene = { i[1]:i[2] for i in exptf}

    ## generate TRN of genes
    trn = []
    for row in net:
        tf = tf_gene[ row[0]]
        tg = tus[ row[1]]
        trn.append( [ tf, tg])
    # remove duplicates
    trn = basic.dedup(trn)
    # remove those with no tf name
    trn = [ i for i in trn if i[0] != '']
    # make a nested list of tf and leader genes
    tfld = [ [ i[0], i[1][0]] for i in trn]
    # remove duplicates
    tfld = basic.dedup(tfld)
    # add new rows per gene for tfs with > 1 gene
    edata = []
    for row in tfld:
        if ',' in row[0]:
            comp = row[0].split(', ')
            for c in comp:
                entry = [ c, row[1]]
                edata.append( entry)
        else:
            entry = row
            edata.append(entry)
    # remove duplicate rows
    dedata = basic.dedup(edata)
    # TFs
    t1 = set( i[0] for i in dedata)
    # leader genes
    t2 = set( i[1] for i in dedata if i[1] not in t1)
    # other targets
    t3 = set( j for i in trn for j in i[1]\
        if j not in (t1|t2))

    # Only keep genes known to have no additional sites
    fun = lambda x: x in genes_wo_additionalsites
    out = {'tf':t1,'leader':t2,'rest':t3}
    out = { k:list(filter(fun,v)) for k,v in out.items()}
    return out

## Function to get E coli's TRN
# in which TFs have been removed from TUs
def get_trn():
    ## regulonDB files
    fnet = os.path.join( path, 'network_tf_tu.txt')
    ftu = os.path.join( path, 'TUSet.txt')
    ftf = os.path.join( path, 'tf_exp_regulon.txt')

    # load network
    net = basic.readcsv( fnet, sep='\t')
    # load info on TUs 
    tus = basic.readcsv( ftu, sep='\t')
    # load info on TFs
    tf_data = basic.readcsv( ftf, sep='\t') 

    # no. of rows in the network
    R = len(net)

    # extract tu names from the network
    tu_names = [ re.sub('\[.+\]','',i[1]) for i in net]

    # get targets' gene names corresponding to above tu names
    target_genes = [ 'NA' for _ in range(R)]
    for x,name in enumerate(tu_names):
        for row in tus:
            if name == row[1]:
                target_genes[x] = row[3]
                break
	
    # get gene names corresponding to the TF's present in the network
    tf_genes = ['NA' for _ in range(R)]
    for x in range(R):
        for row in tf_data:
            if net[x][0] == row[1]:
                tf_genes[x] = row[2]
                break
	
    # bring tf and targets together
    regmap = [ [ tf_genes[r], target_genes[r]] for r in range(R)]
	
    # remove entries with missing tf's gene names or segmented genes 
    regmap = [ i for i in regmap if i[0] != '' and \
	        not any('_' in j for j in i) ]
	
    # remove duplicates & split by tfs
    regmap = { k:[ i[1].split(',') for i in v] for k,v in \
            basic.split_data( basic.dedup(regmap), 0).items()}

    # remove tus which are contained in another tu
    regmap = { k:basic.rm_contained(v) for k,v in regmap.items()}
    
    # remove genes absent from the list of genes without additional sites
    regmap = { k:list( filter( None, [ [ j for j in i if j in \
            genes_wo_additionalsites] for i in v ])) \
            for k,v in regmap.items() if k in\
            genes_wo_additionalsites}
    
    # remove tfs with no genes left
    regmap = { k:v for k,v in regmap.items() if len(v) > 0}
    

    ## remove TFs from TUs
    # TFs
    tfs = regmap.keys()
    # no. of TFs
    nf = len(regmap)
    # initialize final output dictionary
    out = {}
    for f in tfs:
        tus = regmap[f]
        # no. of TUs
        ntu = len(tus)
        # initialize list to hold modified TUs
        modus = [ [] for i in range(ntu)] 
        for i in range(ntu):
            # list of genes in the TU
            tgs = tus[i]
            # remove tfs from the above list of target genes
            modus[i] = [ g for g in tgs if g not in tfs]
        # remove empty TUs
        modus = [ u for u in modus if len(u) >= 1]
        # if there is at least 1 TU left, make an entry in the modified TRN
        if len(modus) >= 1:
            out[f] = modus
    return out

## dict of physico-chemical classes of amino acids 
# and associated function to test whether a change is conservative
from Bio.SeqUtils import IUPACData
aa3to1 = IUPACData.protein_letters_3to1
f_in = os.path.join( path, 'aa_classes.csv')
aa_class = { i[0]:list(filter(None,i[1:])) for i in \
            zip(*basic.readcsv(f_in))}
# convert 3 letter codes to 1 letter
aa_class = { k:[aa3to1[i] for i in v] for k,v in aa_class.items()}
# function
def is_conservative(original,mutant,classdict = aa_class):
    o = original
    m = mutant
    # find original's class
    oc = [ k for k,v in classdict.items() if o in v][0]
    # check if mutant is in the original's class
    if m in classdict[oc]: 
        return True
    else:
        return False

## Make a dict of reference genes with a list of their four-fold degenerate sites
from utilities3.specific import n2c,codon_table
csites = codonsitecount()
ct = codon_table[11]
stops = ct.stop_codons
gene_ffd = {}
for k,v in seqrecord.items():
    # generate codon sequence
    cseq = n2c(v)[:-1]
    # skip the gene if there are any interemediate stop codons 
    if any( i in stops for i in cseq): continue
    # 0-based index of codons with ffd sites
    entry = [ x for x,i in enumerate(cseq) if csites[i][2] == 1]
    gene_ffd[k] = entry

