#!/bin/bash
#SBATCH --job-name=narw_mask_repeats
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=25G
#SBATCH --time=00:30:00

bedtools intersect -v -a ~/projects/def-frasiert/RW_WGS/vcf/NARW/merged_narw_all_scaffold.vcf.gz \
	-b ~/scratch/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz \
	-wa -header > ~/scratch/temp/filtervariants/narw_masked.vcf
