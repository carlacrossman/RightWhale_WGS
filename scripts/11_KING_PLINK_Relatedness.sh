#!/bin/bash
#SBATCH --job-name=PLINK_KING
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00

## From inside ~/scratch/relatedness/
## Run PLINK to generate .bed, .fam and .bim files

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered_aug5.vcf.gz --allow-extra-chr --out narw

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_masked_filtered_aug5.vcf.gz --allow-extra-chr --out srw

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/all_on_narw_masked_filtered_aug5.vcf.gz --allow-extra-chr --out all

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/SRW_w_bowhead_masked_filtered.vcf.gz --allow-extra-chr --out srw_w_bowhead

## Edit the .bim file to remove the HiC_scaffold_ prefix from column 1 and rename a few to clearly identify sex chromosome


sed 's,HiC_scaffold_\(.*\),\1,' all.bim > all_edit.bim
sed 's,HiC_scaffold_\(.*\),\1,' narw.bim > narw_edit.bim
sed 's,HiC_scaffold_\(.*\),\1,' srw_w_bowhead.bim > srw_w_bowhead_edit.bim
sed -i 's,21\(\t.*\),1\1,' srw_w_bowhead_edit.bim
sed -i 's,8\(\t.*\),21\1,' srw_w_bowhead_edit.bim
sed 's,HiC_scaffold_\(.*\),\1,' srw.bim > srw_edit.bim
sed -i 's,21\(\t.*\),1\1,' srw_edit.bim
sed -i 's,8\(\t.*\),21\1,' srw_edit.bim

## Run KING

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/narw.bed --fam ~/scratch/relatedness/narw.fam --bim ~/scratch/relatedness/narw_edit.bim --kinship --degree 3 --prefix narw_kinship --sexchr 21

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/all.bed --fam ~/scratch/relatedness/all.fam --bim ~/scratch/relatedness/all_edit.bim --kinship --degree 3 --prefix all_kinship --sexchr 21

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/srw.bed --fam ~/scratch/relatedness/srw.fam --bim ~/scratch/relatedness/srw_edit.bim --kinship --degree 5 --prefix srw_kinship --sexchr 21

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/srw_w_bowhead.bed --fam ~/scratch/relatedness/srw_w_bowhead.fam --bim ~/scratch/relatedness/srw_w_bowhead_edit.bim --kinship --degree 5 --prefix srw_w_bowhead_kinship --sexchr 21
