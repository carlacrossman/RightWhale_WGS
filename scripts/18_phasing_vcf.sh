#!/bin/bash
#SBATCH --job-name=phasing_vcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=12:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

cp HiC_scaffold_${CHR}_all_on_narw_filtered_aug8* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}

#Set up Input Files
UNPHASED_VCF=HiC_scaffold_${CHR}_all_on_narw_filtered_aug8.vcf.gz
UNPHASED_VCF_NOMAS=HiC_scaffold_${CHR}.noMultiAllelicSites.vcf.gz
LOG_ALIGN=HiC_scaffold_${CHR}.alignments
LOG_MAIN=HiC_scaffold_${CHR}.main

PHASED_HAPS=HiC_scaffold_${CHR}.phased.haps.gz
PHASED_SAMPLE=HiC_scaffold_${CHR}.phased.samples
PHASED_VCF=HiC_scaffold_${CHR}.onlyPhased.vcf
NEW_HEAD=HiC_scaffold_${CHR}_newhead.onlyPhased.vcf

LOG_CONVERT=HiC_scaffold_${CHR}.convert
FINAL_VCF=HiC_scaffold_${CHR}.phased.vcf


#Preparation
bcftools view -M 2 -O z $UNPHASED_VCF > $UNPHASED_VCF_NOMAS
shapeit -check --input-vcf $UNPHASED_VCF_NOMAS --output-log $LOG_ALIGN --thread 7

#Main run
shapeit -V $UNPHASED_VCF_NOMAS --output-max $PHASED_HAPS $PHASED_SAMPLE --output-log $LOG_MAIN --thread 7 --force
shapeit -convert --input-haps $PHASED_HAPS $PHASED_SAMPLE --output-vcf $PHASED_VCF --output-log $LOG_CONVERT --thread 7
echo "SHAPEIT DONE"

#new header on phased vcf
bcftools view -h $UNPHASED_VCF | bcftools reheader -h - $PHASED_VCF -o ${NEW_HEAD}
bcftools view -Oz -o $NEW_HEAD.gz $NEW_HEAD
bcftools index ${NEW_HEAD}.gz
echo "REHEADER DONE"


#Extracting sites in the unphased vcf, not present in the phased sites (i.e. those sites which could not be phased)

bcftools isec -n -1 -p ${CHR}/ -c all $UNPHASED_VCF $NEW_HEAD.gz
echo "ISEC DONE"

#Concatenate unphased sites with phased sites and sort the vcf.
bcftools concat $NEW_HEAD ${CHR}/0000.vcf -Oz -o $FINAL_VCF
echo "CONCAT DONE"
bcftools sort -Oz -o HiC_scaffold_${CHR}_sorted.vcf.gz $FINAL_VCF

mv HiC_scaffold_${CHR}* ${SLURM_SUBMIT_DIR}/output/
