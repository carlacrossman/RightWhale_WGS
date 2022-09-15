#!/bin/bash
#SBATCH --job-name=rw_sortsam
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=10
#SBATCH --cpus-per-task=30
#SBATCH --mem=60G
#SBATCH --time=10:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/rw_fastq_list)

java -Xmx40G -jar ${GATK_JAR} SortSam \
  -I ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  -O ~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  --TMP_DIR ~/scratch/temp_mapping/ \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-SortSam.out

