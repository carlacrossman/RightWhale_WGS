##########################################
###---   Trimming Raw Fastq Files   ---###
##########################################
#
# This script trims the fastq files based on quality scores
# and produces four output files: forward and reverse paired reads (FP & RP)
# forward and reverse unpaired reads (FUP & RUP). 
# Unpaired reads will not be used in subsequent analyses. 
# Standard Error is saved to QC/std_out/.
#
# The script also runs fastqc to generate quality plots after trimming.
#
# Finally this script helps clean up after itself by removing the unpaired
# fastq files from trimmomatic.
#
# This script removes reads with: 
# 	- average quality below Q30
# 	- minimum length less than 36bp
#	- trims leading bases with bases less than Q20
#	- trims in a sliding window of 5bp with mean quality of 20
#	- removes illumina adapters provided by reference/adapters.fa 
# 
# Requirements:
#	- trimmomatic (module: mugqic/trimmomatic/0.36)
#	- java (module: mugqic/java/openjdk-jdk1.8.0_72)
#	- fastqc (module: mugqic/fastqc/0.11.5)
#	- file 'fastq_list' that specifies samples by lane to be trimmed 
#	- fastq files in raw_fastq/
#	- This was run as an array job (1-69)
#
# To be run from: main project directory
#
##########################################


SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_list)

java -Xmx6G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE \
    -threads 12 -phred33 \
    raw_fastq/${SAMPLELANE}_R1_001.fastq.gz \
    raw_fastq/${SAMPLELANE}_R2_001.fastq.gz \
    ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
    ~/scratch/temp/$SAMPLELANE-FUP.fq.gz \
    ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz \
    ~/scratch/temp/$SAMPLELANE-RUP.fq.gz \
    ILLUMINACLIP:reference/adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
    AVGQUAL:30 MINLEN:36 2> QC/std_out/$SAMPLELANE.trim.out

fastqc -o QC/trimQC/ ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz
fastqc -o QC/trimQC/ ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz

rm ~/scratch/temp/$SAMPLELANE-*.fq.gz

