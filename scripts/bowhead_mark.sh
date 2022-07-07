#!/bin/bash
#SBATCH --job-name=bowhead_mark_bam
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=30G
#SBATCH --time=08:00:00

BAMNAME="SRR1685383"

java -Xmx25G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  --INPUT ~/scratch/merged_reads/$BAMNAME-merged_min1Mb.bam \
  --OUTPUT ~/projects/def-frasiert/RW_WGS/merged_bam/$BAMNAME-merged_marked.bam \
  --TMP_DIR ~/scratch/temp_marking/ \
  --METRICS_FILE=QC/mergeQC/$BAMNAME-merged_marked.metrics


