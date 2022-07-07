##########################################
###--- Merging .bam and Marking Dups---###
##########################################
#
# This script merges bam files from 3 lanes into 
# a single merged bam file.  
# Standard Error is saved to QC/std_out/.
#
# The script also runs MarkDuplicates to mark (not remove)
# duplicate reads and index these bam files. 
# Metric files are written to QC/mergeQC/
# 
# Requirements:
#	- java (module: mugqic/java/openjdk-jdk1.8.0_72)
#	- GATK4 (module: mugqic/GenomeAnalysisTK/4.1.0.0)
#	- file 'bam_list' that specifies the sample name prefixes to use
#	- .bam files for each lane in ~/scratch/mapped_reads/
#	- This was run as an array job
#
# To be run from: main project directory
#
##########################################

BAMNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bam_list)
SAMPLE=$(grep -oP '(?<=\b).*?(?=_2)' <<< "$BAMNAME")


java -Xmx25G -jar ${GATK_JAR} MergeSamFiles \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L002-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L003-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L004-aln-sorted_min1Mb.bam \
  --OUTPUT ~/scratch/merged_reads/$SAMPLE-merged_min1Mb.bam \
  2> QC/std_out/$SAMPLE-merge.out

java -Xmx25G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  --INPUT ~/scratch/merged_reads/$SAMPLE-merged_min1Mb.bam \
  --OUTPUT merged_bam/$SAMPLE-merged_marked.bam \
  --TMP_DIR ~/scratch/temp_marking/ \
  --METRICS_FILE QC/mergeQC/$SAMPLE-merged_marked.metrics


