#!/bin/bash
#SBATCH --job-name=BCFTOOLS_ROH
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=12:00:00

VCF=~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz
KEEP=../narw
bcftools roh -e ${KEEP} -G 30 -Orz --threads 8 -o ${VCF}_bcftoolsROH.txt.gz ${VCF}
