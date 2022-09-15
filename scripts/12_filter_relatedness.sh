#!/bin/bash
#SBATCH --job-name=filter_relatedness
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1:00:00



#vcftools --remove-indv EGL312-1a --remove-indv SID181803 --remove-indv EGL276-1 --not-chr HiC_scaffold_21 --not-chr HiC_scaffold_111 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered_aug5.vcf.gz --recode --temp ${SLURM_TMPDIR} --stdout | gzip -c > ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz

vcftools --remove-indv Eau019 --not-chr HiC_scaffold_8 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_masked_filtered_aug5.vcf.gz --temp ${SLURM_TMPDIR} --recode --stdout | gzip -c > ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz

#vcftools  --remove-indv EGL312-1a --remove-indv SID181803 --remove-indv EGL276-1 --remove-indv Eau019 --not-chr HiC_scaffold_21 --not-chr HiC_scaffold_111 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/all_on_narw_masked_filtered_aug5.vcf.gz --temp ${SLURM_TMPDIR} --recode --stdout | gzip -c > ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz

vcftools --remove-indv Eau019 --not-chr HiC_scaffold_8 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/SRW_w_bowhead_masked_filtered.vcf.gz --temp ${SLURM_TMPDIR} --recode --stdout | gzip -c > ~/projects/def-frasiert/RW_WGS/vcf/locked/SRW_w_bowhead_unrelated_filtered_aug5.vcf.gz
