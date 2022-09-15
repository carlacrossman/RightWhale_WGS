#!/bin/bash
#SBATCH --job-name=bowhead_mergemark_bam
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-3,5-22%12
#SBATCH --cpus-per-task=20
#SBATCH --mem=30G
#SBATCH --time=06:00:00

BAMNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/bam_list)
SAMPLE=$(grep -oP '(?<=\b).*?(?=_2)' <<< "$BAMNAME")


java -Xmx25G -jar ${GATK_JAR} MergeSamFiles \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L002-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L003-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L004-aln-sorted_min1Mb.bam \
  --OUTPUT ~/scratch/merged_reads/$SAMPLE-merged_min1Mb.bam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLE-merge.out

java -Xmx25G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  --INPUT ~/scratch/merged_reads/$SAMPLE-merged_min1Mb.bam \
  --OUTPUT ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam \
  --TMP_DIR ~/scratch/temp_marking/ \
  --METRICS_FILE=QC/mergeQC/$SAMPLE-merged_marked.metrics


