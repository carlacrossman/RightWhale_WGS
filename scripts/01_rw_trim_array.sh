#!/bin/bash
#SBATCH --account=def-frasiert
#SBATCH --job-name=trim_sample
#SBATCH --array=1-69%9
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=12
#SBATCH --time=04:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/fastq_list)

java -Xmx6G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE \
    -threads 12 -phred33 \
    ~/projects/def-frasiert/RW_WGS/raw_fastq/${SAMPLELANE}_R1_001.fastq.gz \
    ~/projects/def-frasiert/RW_WGS/raw_fastq/${SAMPLELANE}_R2_001.fastq.gz \
    ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
    ~/scratch/temp/$SAMPLELANE-FUP.fq.gz \
    ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz \
    ~/scratch/temp/$SAMPLELANE-RUP.fq.gz \
    ILLUMINACLIP:reference/adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
    AVGQUAL:30 MINLEN:36 2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE.trim.out

fastqc -o ~/projects/def-frasiert/RW_WGS/QC/trimQC/ ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz
fastqc -o ~/projects/def-frasiert/RW_WGS/QC/trimQC/ ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz

rm ~/scratch/temp/$SAMPLELANE-*.fq.gz

