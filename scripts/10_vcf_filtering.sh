#!/bin/bash
#SBATCH --job-name=vcf_filters
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=1:00:00

## I ran these in an interactive shell, so the file names and directories would need to be updated for future use

bedtools intersect -v -a ~/projects/def-frasiert/RW_WGS/vcf/NARW/merged_narw_all_scaffold.vcf.gz -b ~/scratch/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz -wa -header > ~/scratch/temp/filtervariants/narw_masked.vcf

gzip narw_masked.vcf

vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked.vcf.gz --recode --recode-INFO-all --out ~/scratch/temp/filtervariants/narw_masked_formatfilter_aug5 --minDP 10 --minGQ 30

gzip narw_masked_formatfilter_aug5.recode.vcf

bcftools filter -i 'COUNT(FORMAT/GT="mis")<=3 && INFO/MQ>30 && INFO/DP<790' -Oz -o ~/scratch/temp/filtervariants/narw_masked_filtered_aug5.vcf.gz   ~/scratch/temp/filtervariants/narw_masked_formatfilter_aug5.recode.vcf.gz

