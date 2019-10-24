# generate complete TRN of E coli as dicts of TF with nested lists of TUs

# require
import re
from utilities3 import basic

def get_trn():
    ## regulonDB files
    fnet = 'network_tf_tu.txt'
    ftu = 'TUSet.txt'
    ftf = 'tf_exp_regulon.txt'

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
    return regmap
