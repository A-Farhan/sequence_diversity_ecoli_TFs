# script to submit qsub jobs for processing of SNP matrix

## Arguments
# For SGE
queue="all.q"

# path for scripts
SCRIPTPATH=<full path to the repository>
# full path to input directory
dr="$1"
# path to python program
PY=<python>

# files
f1="genes_present.txt"
f2="$dr/genes_pa.tab"
f3="$dr/snps_filtered.tab"
f4="variants.tab"

# sub-directories of individual sequencing runs for the project
sdirs=(`ls -d $dr/*RR*/`)

# no. of runs
L=${#sdirs[@]}

#clean batch system log files

if [ -s $dr/qsub_ids.txt ]; then
    rm $dr/qsub_ids.txt
fi
if [ -s $dr/patab_id.txt ]; then
    rm $dr/patab_id.txt
fi
if [ -s $dr/filter_id.txt ]; then
    rm $dr/filter_id.txt
fi
if [ -s $dr/annot_id.txt ]; then
    rm $dr/annot_id.txt
fi

wdr=${dr/$HOME\/varspots\/}
echo "directory: $wdr"

# script 1: infer genes presence for each run
for(( i=0; i<L; i++ )); do    
    # extract run accession
    run=$(basename ${sdirs[$i]/\//})
    # log file path
    logf="${sdirs[$i]}pa_genes.log"
    # output file name
    outf="${sdirs[$i]}$f1"   
    # list of variables
    var="dr=${sdirs[$i]},SCRIPTPATH=$SCRIPTPATH"
    if [[ ! -s $outf ]]; then
        echo "submitting qsub job to infer genes presence/absence for $run"
        qsub_id=`qsub -N pa_genes -j y -o $logf -v "$var" $SCRIPTPATH/pa_genes.sh` 
        echo -e "$qsub_id" >> $dr/qsub_ids.txt
    fi
done

# script 2: generate gene P/A table 
if [[ ! -s $f2 ]]; then
    echo "submitting qsub job to make gene presence/absence table"
    if [[ -s $dr/qsub_ids.txt ]]; then
        qsub_cat_ids=`cat $dr/qsub_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
        depend="-hold_jid ${qsub_cat_ids}"
        qsub_id=`qsub $depend gene_pa_table.job $wdr`
        if [[ $? -ne 0 ]]; then 
            qsub_id=`qsub gene_pa_table.job $wdr`
        fi
    else
        qsub_id=`qsub gene_pa_table.job $wdr`
    fi
    echo -e "$qsub_id" >> $dr/patab_id.txt
fi

# script 3: extract and filter SNPs
if [[ ! -s $f3 ]]; then
    echo "Submitting qsub job for SNP filtering"
    if [[ -s $dr/patab_id.txt ]]; then
        qsub_cat_ids=`cat $dr/patab_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
        depend="-hold_jid ${qsub_cat_ids}"
        qsub_id=`qsub $depend filter_snp.job $wdr`
    else
        qsub_id=`qsub filter_snp.job $wdr`
    fi
    echo -e "$qsub_id" >> $dr/filter_id.txt
fi

# script 4: break SNP matrix into codon files of individual runs
if [[ `ls $dr/*RR*/$f4 | wc -l` -lt $L ]]; then
    echo "Submitting qsub job for breaking SNP matrix"
    if [[ -s $dr/filter_id.txt ]]; then
        qsub_cat_ids=`cat $dr/filter_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
        depend="-hold_jid ${qsub_cat_ids}"
        qsub_id=`qsub $depend break_snpmat.job $dr`
    else
        qsub_id=`qsub break_snpmat.job $dr`
    fi
    echo -e "$qsub_id" >> $dr/break_id.txt
fi

exit 0
