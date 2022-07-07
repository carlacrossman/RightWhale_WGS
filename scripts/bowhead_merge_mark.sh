#!/bin/bash
#SBATCH --job-name=bowhead_mergemark_bam_SRW
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=30G
#SBATCH --time=06:00:00

BAMNAME="SRR1685383"

java -Xmx25G -jar ${GATK_JAR} MergeSamFiles \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}.L003-aln-sorted_min1Mb_SRW.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}.L004-aln-sorted_min1Mb_SRW.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}.L006-aln-sorted_min1Mb_SRW.bam \
  --OUTPUT ~/scratch/merged_reads/$BAMNAME-merged_min1Mb_SRW.bam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$BAMNAME-merge_SRW.out

java -Xmx25G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  --INPUT ~/scratch/merged_reads/$BAMNAME-merged_min1Mb_SRW.bam \
  --OUTPUT ~/projects/def-frasiert/RW_WGS/merged_bam/$BAMNAME-merged_marked_SRW.bam \
  --TMP_DIR ~/scratch/temp_marking/ \
  --METRICS_FILE ~/projects/def-frasiert/RW_WGS/QC/mergeQC/$BAMNAME-merged_marked_SRW.metrics


