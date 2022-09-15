#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-frasiert
#SBATCH --job-name=ref_alignment
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8

bwa index reference/Eubalaena_glacialis_HiC_min1Mb.fasta

samtools faidx reference/Eubalaena_glacialis_HiC_min1Mb.fasta

samtools dict reference/Eubalaena_glacialis_HiC_min1Mb.fasta -o reference/Eubalaena_glacialis_HiC_min1Mb.dict

