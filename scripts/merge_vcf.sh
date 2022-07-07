#!/bin/bash
#SBATCH --job-name=merge_vcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=4:00:00

vcf-merge srw_HiC_scaffold_1.vcf.gz srw_HiC_scaffold_2.vcf.gz srw_HiC_scaffold_3.vcf.gz srw_HiC_scaffold_4.vcf.gz srw_HiC_scaffold_5.vcf.gz srw_HiC_scaffold_6.vcf.gz srw_HiC_scaffold_7.vcf.gz srw_HiC_scaffold_8.vcf.gz srw_HiC_scaffold_9.vcf.gz srw_HiC_scaffold_10.vcf.gz srw_HiC_scaffold_11.vcf.gz srw_HiC_scaffold_12.vcf.gz srw_HiC_scaffold_13.vcf.gz srw_HiC_scaffold_14.vcf.gz srw_HiC_scaffold_15.vcf.gz srw_HiC_scaffold_16.vcf.gz srw_HiC_scaffold_17.vcf.gz srw_HiC_scaffold_18.vcf.gz srw_HiC_scaffold_19.vcf.gz srw_HiC_scaffold_20.vcf.gz srw_HiC_scaffold_21.vcf.gz | bgzip -c > merged_srw_all_scaffold.vcf.gz


