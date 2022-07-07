#!/bin/bash
#SBATCH --job-name=srw_mask_repeats
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=25G
#SBATCH --time=00:30:00

bedtools intersect -v -a ~/projects/def-frasiert/RW_WGS/vcf/SRW/merged_srw_all_scaffold.vcf.gz \
	-b ~/scratch/RWref_HiC.repeatmasker.trf.windowmasker.gff.gz \
	-wa -header > ~/scratch/temp/filtervariants/srw_masked.vcf
