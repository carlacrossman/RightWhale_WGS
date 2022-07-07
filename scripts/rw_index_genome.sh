##########################################
###---Indexing NARW Reference Genome---###
##########################################
#
# This script creates the required index files
# from the NARW genome for subsequent analyses. 
# 
# Requires: reference genome saved as: reference/Eubalaena_glacialis_HiC_min1Mb.fasta 
#
# To be run from: main project directory
#
##########################################


bwa index reference/Eubalaena_glacialis_HiC_min1Mb.fasta

samtools faidx reference/Eubalaena_glacialis_HiC_min1Mb.fasta

samtools dict reference/Eubalaena_glacialis_HiC_min1Mb.fasta -o reference/Eubalaena_glacialis_HiC_min1Mb.dict

