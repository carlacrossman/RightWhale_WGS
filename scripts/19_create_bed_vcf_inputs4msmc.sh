#!/bin/bash
#SBATCH --job-name=create_bed_vcf_inputs
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-189%10
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=1:00:00

################################
### The purpose of this script is to create the bed files required for generate_multihetsep.py for MSMC
###
### Requirements:
###		1) A file with two columns with each row being unique combination of sample name and scaffold 
###		2) A gvcf file for each individual, at each scaffold
###		3) Python script installed gvcf2bed
###		4) A phased vcf for each chromosome
###		5) A vcf of raw variants created after joint variant calling before filtering
###		6) vcftools/0.1.14
###		7) bcftools/1.11
###		8) bedtools/2.30.0
###
################################


###---INPUT FILES AND VARIABLES---###

# Individual and Scaffold/Chromosome List
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" NARW_IND_SCAFFOLDS)
CHR=$(echo $FILE | awk '{print $1}')
SAMPLE=$(echo $FILE | awk '{print $2}')

# Parent directory for gvcf files (Mine are subset with species and scaffold, adjust below as needed)
GVCF_FILE=~/scratch/gvcf/SRWonNARW/HiC_scaffold_${CHR}/${SAMPLE}-HiC_scaffold_${CHR}.g.vcf.gz

# VCF file after joint variant calling, before filtering variants
RAW_VCF=~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz

# VCF file after all filtering and phasing - the VCF file you deemed 'locked in'
PHASED_VCF=~/scratch/msmc/shapeit/phased/HiC_scaffold_${CHR}_sorted.vcf.gz

# OUTPUT DIRECTORY
OUT_DIR=~/scratch/msmc/generate_multisephet_input/MSMC_IM_runs/

###---RUN CODE IN AN ARRAY WITH EACH JOB A DIFFERENT COMBINATION OF SCAFFOLD AND CHROMOSOME ---###

#create bed file from gvcf file to mark all sites that were able to be mapped
gvcf2bed -I ${GVCF_FILE} -O ${SAMPLE}_${CHR}.bed

#split the phased vcf file into individuals, remove indels and missing data.
vcftools --gzvcf ${PHASED_VCF} --indv ${SAMPLE} --max-missing 1.0 --remove-indels --recode --out ${SAMPLE}_${CHR}_all.filtered

#split the raw vcf file by individual and chromosome
vcftools --gzvcf ${RAW_VCF} --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}_all.originalmerge

#create a vcf file of sites where we called variants, but that were filtered out for any reason
bedtools subtract -a ${SAMPLE}_${CHR}_all.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}_all.filtered.recode.vcf > ${SAMPLE}_${CHR}_all_filtered_sites.vcf

#correct header
bcftools view -h ${SAMPLE}_${CHR}_all.filtered.recode.vcf > temp-${SAMPLE}-${CHR}_all
cat temp-${SAMPLE}-${CHR}_all ${SAMPLE}_${CHR}_all_filtered_sites.vcf > ${SAMPLE}_${CHR}_all_filtered_sites_wheader.vcf

#use bedtools to create a bed file that removed the sites that were filtered out from the original bed file that listed all sites we could map as we are not confident these sites are not variants and we don't want to make this assumption
bedtools subtract -a ${SAMPLE}_${CHR}.bed -b ${SAMPLE}_${CHR}_all_filtered_sites_wheader.vcf > ${SAMPLE}_${CHR}_all_msmcmask.bed
gzip ${SAMPLE}_${CHR}_all_msmcmask.bed

#remove vcf file with raw variant information as it is no longer needed
rm ${SAMPLE}_${CHR}_all.originalmerge.recode.vcf
#remove tempory header and temporary vcf without correct header
rm temp-${SAMPLE}-${CHR}_all
rm ${SAMPLE}_${CHR}_all_filtered_sites.vcf

#Move good output files
mv ${SAMPLE}_${CHR}_all.filtered.recode.vcf ${OUT_DIR}
mv ${SAMPLE}_${CHR}_all_msmcmask.bed.gz ${OUT_DIR}
