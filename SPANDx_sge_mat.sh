#!/bin/bash

# only run function: Matrix and MergeBED, assuming variant calling is over

#####################################################
# 
# Thanks for using SPANDx!!
#
# USAGE: SPANDx.sh <parameters, required> -r <reference, without .fasta extension> [parameters, optional] -o [organism] -m [generate SNP matrix yes/no] -i [generate indel matrix yes/no] 
# -a [include annotation yes/no] -v [Variant genome file. Name must match the SnpEff database] -s [Specify read prefix to run single strain or none to construct a SNP matrix from a previous analysis ] -t [Sequencing technology used Illumina/Illumina_old/454/PGM] 
# -p [Pairing of reads PE/SE] -w [Window size in base pairs for BEDcoverage module]
#
# SPANDx by default expects reads to be paired end, in the following format: STRAIN_1_sequence.fastq.gz for the first pair and STRAIN_2_sequence.fastq.gz for the second pair. 
# Reads not in this format will be ignored.
# If your data are not paired, you must set the -p parameter to SE to denote unpaired reads. By default -p is set to PE.
#
# SPANDx expects at least a reference file in FASTA format. 
# For compatibility with all programs in SPANDx, FASTA files should conform to the specifications listed here: http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml
# Note that the use of nucleotides other than A, C, G, or T is not supported by certain programs in SPANDx and these should not be used in reference FASTA files. 
# In addition, Picard, GATK and SAMtools handle spaces within contig names differently. Therefore, please avoid using spaces or special characters (e.g. $/*) in contig names.
# 
# By default all read files present in the current working directory will be processed. 
# Sequence files within current directory will be aligned against the reference using BWA, SNPs and indels will be called with GATK and a SNP 
# matrix will be generated with GATK and VCFTools
# 
# Optionally to run the SNP matrix without processing any sequence files set -s to none and -m to yes. If -s is set to none SPANDx will skip all sequence files in the current directory and will not perform the alignment or variant calls. Instead SPANDx will merge all VCF files contained in $PWD/Phylo/snps, interrogate those calls in the alignment files contained within $PWD/Phylo/bams and output a SNP matrix in $PWD/Phylo/out. Before running this module please check that all VCF files located within $PWD/Phylo/snps match the alignment files within $PWD/Phylo/bams
# e.g. SPANDx.sh -r K96243 -m yes -s none
#
# Written by Derek Sarovich and Erin Price - Menzies School of Health Research, Darwin, Australia
# Please send bug reports to mshr.bioinformatics@gmail.com
# If you find SPANDx useful and use it in published work please cite - SPANDx: a genomics pipeline for comparative analysis of large haploid whole genome re-sequencing datasets - BMC Research Notes 2014, 7:618"
# Version 2.2
# 2.0-2.1 Added SGE job handling
# 2.1-2.2 Added SLURM job handling
# 2.3 Added "no resource manager (NONE)" for job handling
#
#################################################################

usage()
{
echo -e  "USAGE: SPANDx.sh <parameters, required> -r <reference, without .fasta extension> [parameters, optional] -o [organism] -m [generate SNP matrix yes/no] -i [generate indel matrix yes/no] -a [include annotation yes/no] -v [Variant genome file. Name must match the SnpEff database] -s [Specify read prefix to run single strain or none to just construct SNP matrix] -t [Sequencing technology used Illumina/Illumina_old/454/PGM] -p [Pairing of reads PE/SE] -w [BEDTools window size in base pairs] -z [include tri-allelic and tetra-allelic SNPs yes/no]"
}
help()
{
usage
cat << _EOF_ 

Thanks for using SPANDx!!

SPANDx requires a reference in FASTA format and for this to be specified with the -r flag. Please do not include the .fasta extension in the reference name.

FASTA files should conform to the specifications listed here: http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml. Some programs used within SPANDx do not allow the use of IUPAC codes and these should not be used in reference FASTA files.

To output a merged SNP file for phylogenetic reconstruction, set the -m flag to yes.

By default, SPANDx will process all sequence data in the present working directory.

By default, SPANDx expects reads to be paired-end (PE) Illumina data, named in the following format: STRAIN_1_sequence.fastq.gz for the first pair and STRAIN_2_sequence.fastq.gz for the second pair. Reads not in this format will be ignored by SPANDx.

If you do not have PE data, you must set the -p flag to SE to denote single-end reads. By default, -p is set to PE.

Sequences will be aligned against the reference using BWA. SNPs and indels will be called with GATK and a SNP matrix will be generated with GATK and VCFTools.

For a more detailed description of how to use SPANDx and its capability, refer to the SPANDx manual or the BMC Research Notes publication. If you have any questions, requests or issues with SPANDx, please send an e-mail to mshr.bioinformatics@gmail.com or derek.sarovich@menzies.edu.au.

If you use SPANDx in published work please cite it!: SPANDx: a genomics pipeline for comparative analysis of large haploid whole genome re-sequencing datasets - BMC Research Notes 2014, 7:618

_EOF_

}

#Define path to SPANDx install
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

# source dependencies
source "$SCRIPTPATH"/SPANDx.config 
source "$SCRIPTPATH"/scheduler.config

#declare variables
declare -rx SCRIPT=${0##*/}

OPTSTRING="o:d:m:a:v:s:t:p:w:i:z:"

declare SWITCH
org=haploid
seq_directory="$PWD"
matrix=yes
annotate=no
strain=all
tech=Illumina
pairing=PE
variant_genome=no
window=1000
indel_merge=no
tri_tetra_alelic=no

# Examine individual options
while getopts "$OPTSTRING" SWITCH; do 
		case $SWITCH in
		       
		o) org="$OPTARG"
		   echo "Organism = $org"
		   ;;
		
		d) seq_directory="$OPTARG"
		   echo "SPANDx expects read files to be in $seq_directory"
		   ;;
		   
        m) matrix="$OPTARG"
           if [ "$matrix" == yes -o "$matrix" == no ]; then 
           echo "SNP matrix will be generated? $matrix"
		   else
		       echo -e "\nIf the optional -m parameter is used it must be set to yes or no\n"
			   echo -e "By default -m is set to no and SPANDx will not generate a SNP matrix\n"
		       usage
			   exit 1
           fi		
           ;;  
		   
		a) annotate="$OPTARG"   
		   if [ "$annotate" == yes -o "$annotate" == no ]; then
		   echo "VCF files will be annotated with SnpEff? $annotate"
		   else
		       echo -e "\nIf the optional -a parameter is used it must be set to yes or no\n"
			   echo -e "By default -a is set to no and variants will not be annotated\n\n"
		       usage
			   exit 1
           fi		
           ;; 
		
		s) strain="$OPTARG"
		   echo "Only strain $strain will be processed. By default all strains within current directory will be processed"
		   ;;
		
		t) tech="$OPTARG"
			if [ "$tech" == Illumina -o "$tech" == Illumina_old -o "$tech" == 454 -o "$tech" == PGM ]; then
			    echo -e "Technology used is $tech\n"
			else
			   echo -e "If the optional parameter -t is used it must be set to Illumina, Illumina_old, 454 or PGM.\n"
			   echo -e "By default -t is set to Illumina and SPANDx will assume your sequence data is in Illumina format with quality encoding of 1.8+\n"
			   usage
			   exit 1
			fi
			;;
			
		p) pairing="$OPTARG"
		   if [ "$pairing" == PE -o "$pairing" == SE ]; then 
               echo "The pairing of your reads has been indicated as $pairing"
		   else
		       echo -e "\nIf the optional -p parameter is used it must be set to PE or SE\n"
			   echo -e "By default -p is set to PE and SPANDx assumes your reads are paired end\n"
		       usage
			   exit 1
		   fi
		   ;;
		   
		w) window="$OPTARG"
		   echo "BEDTools BEDCoverage window size has been set to $window base pairs. Default is 1000bp"
		   ;;
		
		v) variant_genome="$OPTARG"
		   echo "SnpEff will use $variant_genome as the annotated reference"
		   echo "Please make sure this name matches the name in SnpEff database. Also verify chromosome names are correct. These can be found easily at ftp://ftp.ncbi.nih.gov/genomes"
           ;;
		
        i) indel_merge="$OPTARG"
		   if [ "$indel_merge" == yes -o "$indel_merge" == no ]; then
               echo -e "Indels will be merged and checked across genomes = $indel_merge\n"  
		   else
		       echo -e "Indel merge must be set to yes or no. Please refer to the manual for more details\n"
			   exit 1
			   fi
           ;;
		   
		z) tri_tetra_allelic="$OPTARG"
           if [ "$tri_tetra_allelic" == yes -o "$tri_tetra_allelic" == no ]; then
               echo -e "Tri- and tetra-allelic SNPs will be included in the phylogenetic analyses\n"  
		   else
		       echo -e "Tri- and tetra-allelic SNPs (-z) must be set to yes or no. Please refer to the manual for more details\n"
			   exit 1
		   fi
           ;;

		
		\?) usage
		    exit 1
		    ;;
		
		h) help
		   exit 1
		   ;;
		   
		*) echo "script error: unhandled argument"
           usage
		   exit 1
		   ;;
		   
		
	esac
done
  
echo -e "\nThe following parameters will be used\n"
echo -e "-------------------------------------\n"
echo "Organism = $org"
echo "Output directory and directory containing sequence files = $seq_directory"
echo "SNP matrix will be created? = $matrix"
echo "Genomes will be annotated? = $annotate"
echo "Strain(s) to be processed = $strain"
echo "Sequence technology used = $tech"
echo "Pairing of reads = $pairing"
echo "Variant genome specified for SnpEff = $variant_genome"
echo "Window size for BedTools = $window"
echo "Indels will be merged and corrected = $indel_merge"
echo -e "-------------------------------------\n\n"

ref='NC_000913'
ref_fasta="$seq_directory/$ref.fasta" # reference file in fasta format
ref_index=$seq_directory/$ref.bwt #index file created with BWA
REF_INDEX_FILE=$seq_directory/$ref.fasta.fai #index created with SAMtools
REF_BED=$seq_directory/$ref.bed
REF_DICT=$seq_directory/$ref.dict #Dictionary file created with Picard tools

if [ ! $PBS_O_WORKDIR ]; then
        PBS_O_WORKDIR="$seq_directory"
fi

cd $PBS_O_WORKDIR

#to do. Needs to have an n > 1 check before the test	
	
if [ "$matrix" == no -a "$indel_merge" == yes ]; then
    echo -e "Indel merge has been requested. As this cannot be determine without creation of a SNP matrix the matrix variable has been set to yes\n"
	matrix=yes
fi
	
## create directory structure
	
if [ ! -d "BEDcov" ]; then
    mkdir $seq_directory/BEDcov 
fi
if [ ! -d "tmp" ]; then
	mkdir $seq_directory/tmp
fi
if [ ! -d "Outputs" ]; then
	mkdir $seq_directory/Outputs
fi
if [ ! -d "Outputs/SNPs_indels_PASS" ]; then
	mkdir $seq_directory/Outputs/SNPs_indels_PASS
fi
if [ ! -d "Outputs/SNPs_indels_FAIL" ]; then
	mkdir $seq_directory/Outputs/SNPs_indels_FAIL
fi
if [ ! -d "Outputs/Comparative" ]; then
	mkdir $seq_directory/Outputs/Comparative
fi
if [ ! -d "logs" ]; then
	mkdir $seq_directory/logs
fi

## checking variables for the annotation module

if [ "$annotate" == yes ]; then
    grep "$variant_genome" "$SNPEFF_CONFIG" &> /dev/null
    status=$?
    if [ ! $status == 0 ]; then
        echo "SPANDx couldn't find the reference genome in the SnpEff config file" 
		echo "The name of the annotated reference genome specified with the -v switch must match a reference genome in the SnpEff database"
        echo "Does the SnpEff.config file contain the reference specified with the -v switch?"
		echo "Is the SnpEff.config file in the location specified by SPANDx.config?"
	    echo "If both of these parameters are correct please refer to the SnpEff manual for further details on the setup of SnpEff"
		exit 1
    else
        echo -e "SPANDx found the reference file in the SnpEff.config file\n" 
    fi
	
	#test to see if the chromosome names in the SnpEff database match those of the reference file
	
	CHR_NAME=`$JAVA -jar $SNPEFF dump "$variant_genome" | grep -A1 'Chromosomes names' | tail -n1 | awk '{print $2}'|sed "s/'//g"`
	REF_CHR=`head -n1 $ref_fasta | sed 's/>//'`  
	if [ "$CHR_NAME" == "$REF_CHR" ]; then
	    echo -e "Chromosome names in the SnpEff database match the reference chromosome names, good\n"
	else
	    echo -e "Chromosome names in the SnpEff database DON'T match the reference chromosome names.\n"
		echo -e "Please change the names of the reference file to match those in the SnpEff database.\n"
		echo -e "If you are unsure what these are, run: $JAVA -jar $SNPEFF dump $variant_genome\n"
		echo -e "The first chromosome name is $CHR_NAME.\n\n"
        echo -e "The reference chromosome name is $REF_CHR.\n\n"
                echo -e "If you choose to continue the annotation component of SPANDx may fail.\n"
                echo -e "Do you want to continue?\n"
                echo -e "Type Y to continue or anything else to exit\n"
                read ref_test
                if [ "$ref_test" == "Y" -o "$ref_test" == "y" -o "$ref_test" == "yes" ]; then
                   echo -e "Continuing\n\n"
                else
                   exit 1
                fi
	fi	
	
	
	if [ ! -d "$SNPEFF_DATA/$variant_genome" ]; then
	    echo -e "Downloading reference genome to SnpEff database\n"
		echo -e "If the program hangs here please check that the proxy settings are correct and the cluster has internet access\n"
		echo -e "If required SNPEff databases can be manually downloaded and added to the SPANDx pipeline\n"
		echo -e "Running the following command:"
		echo "$JAVA $JAVA_PROXY -jar $SNPEFF download -v $variant_genome"
        echo -e "In the following directory $PBS_O_WORKDIR\n"		
		$JAVA ${JAVA_PROXY} -jar $SNPEFF download -v $variant_genome
	else 
        echo -e "Annotated reference database has already been downloaded for SnpEff\n"
    fi	
fi


if [ "$SCHEDULER" == PBS -o "$SCHEDULER" == SGE -o "$SCHEDULER" == SLURM -o "$SCHEDULER" == NONE ]; then
	echo -e "SPANDx will use $SCHEDULER for resource management\n"
	else
	echo -e "SPANDx requires you to set the SCHEDULER variable to one of the following: PBS, SGE, SLURM or NONE. \nIt looks like you might have specified the variable incorrectly.\n"
fi

#clean batch system log files
if [ -s qsub_ids.txt ]; then
    rm qsub_ids.txt
fi
if [ -s qsub_ids2.txt ]; then
    rm qsub_ids2.txt
fi
if [ -s qsub_array_ids.txt ]; then
    rm qsub_array_ids.txt
fi
if [ -s mastervcf_id.txt ]; then
    rm mastervcf_id.txt
fi
if [ -s clean_vcf_id.txt ]; then
	rm clean_vcf_id.txt
fi
if [ -s matrix_id.txt ]; then
    rm matrix_id.txt
fi

# array of sequences to process
sequences_tmp=(`find Outputs/SNPs_indels_PASS/*.snps.PASS.vcf -printf "%f "`)
sequences=("${sequences_tmp[@]/.snps.PASS.vcf/}")

# Functions:
## Matrix is the main comparative genomics section of SPANDx that relies on the outputs from the variant function above
matrix_sge ()			
{	
if [ -s qsub_array_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    qsub_cat_ids=`cat qsub_array_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
    echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/Master_vcf.sh`
	echo -e "$qsub_matrix_id" >> mastervcf_id.txt	
fi
if [ ! -s qsub_array_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/Master_vcf.sh`
	echo -e "$qsub_matrix_id" >> mastervcf_id.txt
fi

### creates clean.vcf files for SNP calls across all genomes

# if a master vcf job has been submitted
if [ -s mastervcf_id.txt ]; then
    # extract its id
    qsub_cat_ids=`cat mastervcf_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	# use it to schedule clean vcf job
    depend="-hold_jid ${qsub_cat_ids}"
	echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	# generate names of output files from vcf file names
    out=("${sequences_tmp[@]/.snps.PASS.vcf/.clean.vcf}")
	# generate bam file names from sequence file names
    bam=("${sequences_tmp[@]/.snps.PASS.vcf/.bam}")
	# store all bam file names in an array after changing their paths to Phylo bam directory
    bam_array=("${bam[@]/#/$PBS_O_WORKDIR/Phylo/bams/}")
    # length of array
    n=${#bam_array[@]}
    # for each bam file stored in the above array
    for (( i=0; i<n; i++ )); do
        # if clean vcf is not present
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${sequences[$i]}.clean.vcf ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/$ref.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
        # Above If indel merging is requested
		if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${sequences[$i]}.clean.vcf -a "$indel_merge" == yes ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi

# Above if master vcf job was not submitted
if [ ! -s mastervcf_id.txt ]; then
    echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	out=("${sequences_tmp[@]/.snps.PASS.vcf/.clean.vcf}")
	bam=("${sequences_tmp[@]/.snps.PASS.vcf/.bam}")
	bam_array=("${bam[@]/#/$PBS_O_WORKDIR/Phylo/bams/}")
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${sequences[$i]}.clean.vcf ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/$ref.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
		if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${sequences[$i]}.clean.vcf -a "$indel_merge" == yes ]; then
		cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/$ref.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
		qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
		echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done 
fi

## if indels merge is set to yes creates clean.vcf files for all indels identified across all genomes

if [ -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    qsub_cat_ids=`cat clean_vcf_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge,tri_tetra_allelic=$tri_tetra_allelic"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$qsub_matrix_id" >> matrix_id.txt
fi
if [ ! -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge,tri_tetra_allelic=$tri_tetra_allelic"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$qsub_matrix_id" >> matrix_id.txt
fi
}

## This function takes the output of each bedcoverage assessment of the sequence alignment and merges them in a comparative matrix
## This function is run when $strain=all but not when a single strain is analysed i.e. $strain doesn't equal all
merge_BED_sge ()
{
if [ -s qsub_array_ids.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt ]; then
    qsub_cat_ids=`cat qsub_array_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
    echo -e "Submitting qsub job for BEDcov merge\n"
    var="seq_path=$seq_directory"
	qsub_BED_id=`qsub -N BEDcov_merge -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/BedCov_merge.sh`
	#echo -e "$qsub_BED_id" >> qsub_BED_id.txt
fi
if [ ! -s qsub_array_ids.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt ]; then
    echo -e "Submitting qsub job for BEDcov merge\n"
    var="ref=$ref,seq_path=$seq_directory"
	qsub_BED_id=`qsub -N BEDcov_merge -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/BedCov_merge.sh`
	#echo -e "$qsub_BED_id" >> qsub_BED_id.txt
fi
}

# run above functions
matrix_sge
merge_BED_sge

exit 0
