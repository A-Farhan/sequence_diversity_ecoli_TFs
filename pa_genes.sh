# infer genes presence/absence based on depth of coverage of each gene region

# when arguments are passed directly to script
dr=$1
SCRIPTPATH=<path to the repository>

# path for programs
Bt=<bedtools_path>
St=<samtools_path>
Rp=<Rscript_path>

# extract sequencing run accession 
run=`basename ${dr%/}`

# get breadth of coverage
INF="$dr/unique/${run}.realigned.bam"
OUTF="$dr/cov_breadth.bed"

if [[ ! -s $OUTF ]]; then
    echo "running command to get breadth of coverage for gene regions"
    $Bt coverage -bed -a $SCRIPTPATH/genes.bed -b $INF > $OUTF
fi

# get depth of coverage
OUTF="$dr/cov_depth.bed"
if [[ ! -s $OUTF ]]; then
    echo "running command to get depth of coverage for gene regions"
    $St bedcov $SCRIPTPATH/genes.bed $INF -Q 50 > $OUTF
fi

# extract lists of genes
OUTF="$dr/genes_present.txt"
if [[ ! -s $OUTF ]]; then
    echo "running R script to make lists of genes present, absent and genes with excess coverage" 
    $Rp $SCRIPTPATH/infer_gene_presence.R $dr
fi
