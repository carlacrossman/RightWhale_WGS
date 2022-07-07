###########################################
###---  Creating genomics databases  ---###
###########################################
#
# This script uses GATK4 GenomicsDBImport to create
# a database of variants across individuals. 
# The path to the genomicsdb-workspace-path argument 
# must be to a non-existant directory. The program manual 
# suggests this can be an empty directory, but it 
# performed better when it created a new directory.
#
# Unlike most scripts, this couldn't be run from the main project directory.
# The paths in the sample-name-map had to be relative paths.
# The genomicsDBimport scripts need to be run from the directory with gvcf files.
# Be sure to update the absolute paths as needed!
# 
# Requirements:
#	- java (module: mugqic/java/openjdk-jdk1.8.0_72)
#	- GATK4 (module: mugqic/GenomeAnalysisTK/4.1.0.0)
#	- file 'HiC_scaffold_19_SRW_samples' that specifies files and samples
#	- This was run as an array job with a job for each scaffold
#
# To be run from: ~/scratch/gvcf/SRW/HiC_scaffold_19/
#
##########################################

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenomicsDBImport \
    --genomicsdb-workspace-path ~/scratch/genomicsDB/SRW_workspace/scaff19/ \
    -L HiC_scaffold_19 \
    --sample-name-map ~/scratch/gvcf/SRW/HiC_scaffold_19/HiC_scaffold_19_SRW_samples \
    --tmp-dir ~/scratch/temp \
    --reader-threads 5
