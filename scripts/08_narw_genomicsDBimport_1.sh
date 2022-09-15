#!/bin/bash
#SBATCH --job-name=narw_genomicsdb_1
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=10:00:00

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenomicsDBImport \
    --genomicsdb-workspace-path ~/scratch/genomicsDB/NARW_workspace/scaff1/ \
    -L HiC_scaffold_1 \
    --sample-name-map ~/scratch/gvcf/NARW/HiC_scaffold_1/HiC_scaffold_1_NARW_samples \
    --tmp-dir ~/scratch/temp \
    --reader-threads 5

