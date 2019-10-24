# Select Data for 
# organism : E coli
# Sequencing platform: Illumina
# library layout: Paired
# library strategy: WGS
# library source: Genomics
# library selection: Random | unspecified | size fractionation | random pcr | pcr

# And Coverage. >= 50X for isolates 

# requires
import os, sys, re
from utilities3.basic import *

# files and paths and variables
fl1 = sys.argv[1]
fl2 = re.sub('.txt','_filtered.txt',fl1)
fl3 = re.sub('.txt','_adrs.txt',fl1)

Cov = int(sys.argv[2]) # coverage cut-off
L = 4641652 # reference genome length for K12 MG1655 NC_000913.3

print("Coverage cut-off set to %i\n"%Cov)

# load data
data1 = readcsv(fl1,'\t')
header = data1[0]
data1 = data1[1:]

print("Total samples: %i\n"%len(data1))

cols = ['scientific_name','instrument_platform','library_layout',\
        'library_strategy','library_selection',\
        'base_count','fastq_ftp']
cxs = [ x for x, i in enumerate(header) if i in cols ]

# filter data by above criteria
data2 = [ i for i in data1 if "Escherichia coli" in i[cxs[0]] and\
        ''.join(i[x] for x in cxs[1:4])=="ILLUMINAPAIREDWGS" \
        and i[cxs[4]] in ['RANDOM','unspecified','size fractionation','RANDOM PCR','PCR'] ]

print("Relevant samples: %i\n"%len(data2))

# filter by coverage
# first, remove missing basecount entries

# get ftp addresses after applying coverage cut-off
ftps,covs,ixs = [],[],[]
for x,row in enumerate(data2):
    bases = row[cxs[5]]
    adrs = row[cxs[6]].split(';')
    # if base count is missing, skip the record
    if bases == '': continue
    # if > 1 ftp addresses are present, then keep only those
    # with end tags
    if len(adrs)>1:
        adrs = [ i for i in adrs if '_' in i ]
    # if coverage satisfies cut-off, take ftp address
    cov = int(bases)/L
    covs.append(cov)
    if cov > Cov: 
        ftps.extend(adrs)
        ixs.append(x)

print("Samples satisfying coverage cut-off: %i\n"%len(ixs))

# write addresses to file
with open(fl3,'w') as flob:
    for i in ftps: flob.write(i+'\n')

# write data left after filtering to another file
out = [ data2[x] for x in ixs ]
writecsv(fl=fl2,data=out,header=header,sep='\t')
