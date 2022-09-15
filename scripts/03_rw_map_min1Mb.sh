#!/bin/bash
#SBATCH --job-name=rw_map_reads
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=28-30,61-69
#SBATCH --cpus-per-task=28
#SBATCH --mem=15G
#SBATCH --time=04:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/fastq_list)
SAMPLE=$(grep -oP '(?<=\b).*?(?=_2)' <<< "$SAMPLELANE")
LANE=${SAMPLELANE: -4}
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:RW\tPU:${SAMPLE}\tCN:MCGILL\tPL:ILLUMINA"

bwa mem -M -t 28 \
  -R $READGROUP \
  ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz > \
  ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-bwa.err

