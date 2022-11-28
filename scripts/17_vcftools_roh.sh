#!/bin/bash
#SBATCH --job-name=VCFTOOLS_ROH
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --mem=10G
#SBATCH --time=12:00:00

VCF=~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz
for CHR in HiC_scaffold_{1..22} ; do
vcftools --gzvcf ${VCF} --LROH --chr ${CHR} --out narw_vcftools_roh_${CHR}
done

