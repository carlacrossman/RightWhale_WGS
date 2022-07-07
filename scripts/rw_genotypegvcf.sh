############################################
###---Genotyping gVCF files by scaffold---###
############################################
#
# This script uses GATK4 GenotypeGVCFs
# to call haplotypes using the genomicsDB created from across all individuals
# Standard Error is saved to QC/gvcfQC/.
#
# Unlike most scripts, this couldn't be run from the main project directory.
# The path to the databses (gendb://$FOLDER) needs to be a relative path.
# Be sure to update the absolute paths as needed!
# 
# Requirements:
#	- java (module: mugqic/java/openjdk-jdk1.8.0_72)
#	- GATK4 (module: mugqic/GenomeAnalysisTK/4.1.0.0)
#	- file 'narw_genDB_inputs' that specifies folders and scaffolds
#	- reference genome .fasta
#	- genomics database for each scaffold indicated by gendb://$FOLDER with a relative path
#	- This was run as an array job with a job for each scaffold
#
# To be run from: ~/scratch/genomicsDB/NARW_workspace/
#
##########################################

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" narw_genDB_inputs)

FOLDER=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenotypeGVCFs \
   -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
   -V gendb://$FOLDER/ \
   -O ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_HiC_scaffold_$SCAFFOLD.vcf.gz \
   2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/narw_HiC_scaffold_$SCAFFOLD-gvcf.out
