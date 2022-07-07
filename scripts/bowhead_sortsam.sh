##########################################
###---Sorting NARW .sam file to .bam---###
##########################################
#
# This script sorts the sam file and saves it
# as a bam files and generates and index.  
# Standard Error is saved to QC/std_out/.
# 
# This code uses the bowhead_fastq_list input file to extract 
# sample name and lane number to use in mapping.
#
# The .sam file should be erased after this step has been completed.
# 
# Requirements:
#	- java (module: mugqic/java/openjdk-jdk1.8.0_72)
#	- GATK4 (module: mugqic/GenomeAnalysisTK/4.1.0.0)
#	- This job required a lot of memory. I used 60G/core
#	- file 'bowhead_fastq_list' that specifies samples by lane to be sorted
#	- mapped .sam files in ~/scratch/temp_mapping/
#
# To be run from: main project directory
#
##########################################

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bowhead_fastq_list)

java -Xmx40G -jar ${GATK_JAR} SortSam \
  -I ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  -O ~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  --TMP_DIR ~/scratch/temp_mapping/ \
  2> QC/std_out/$SAMPLELANE-SortSam.out
