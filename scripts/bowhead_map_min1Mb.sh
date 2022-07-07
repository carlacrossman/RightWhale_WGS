#############################################
###---Mapping bowhead to NARW reference---###
#############################################
#
# This script maps the bowhead fastq reads to the
# scaffolds of the Eubalaena glacialis reference  
# genome that are 1Mbp of longer.  It writes a 
# .sam file to sort in the next step.
# Standard Error is saved to QC/std_out/.
# 
# This code uses the fastq_list input file to extract 
# sample name and lane number to use in mapping.
# 
# Requirements:
#	- bwa (module: mugqic/bwa/0.7.17)
#	- file 'bowhead_fastq_list' that specifies samples by lane to be mapped 
#	- trimmed reads in ~/scratch/trimmed_reds/
#	- This was run as an array job
#
# To be run from: main project directory
#
##########################################

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bowhead_fastq_list)
SAMPLE="SRR1685383"
LANE=${SAMPLELANE: -4}B
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:bowhead\tCN:BOWHEADGENOMEPROJECT\tPL:ILLUMINA"

bwa mem -M -t 28 \
  -R $READGROUP \
  reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz > \
  ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  2> QC/std_out/$SAMPLELANE-bwa.err