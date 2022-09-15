#!/bin/bash
#SBATCH --job-name=qualimap_testing
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=30
#SBATCH --mem=30G
#SBATCH --time=48:00:00

qualimap multi-bamqc -d qualimap_sample_file \
    -outdir ~/projects/def-frasiert/QC/mapQC -nr 10000 \
    --java-mem-size=25G -r
