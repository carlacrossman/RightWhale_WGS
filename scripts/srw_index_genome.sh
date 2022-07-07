##########################################
###--- Indexing SRW Reference Genome---###
##########################################
#
# This script creates the required index files
# from the SRW genome for subsequent analyses. 
# 
# Requires: reference genome saved as: reference/RWref_HiC_min1Mb.fasta
#
# To be run from: main project directory
#
##########################################

bwa index reference/RWref_HiC_min1Mb.fasta

samtools faidx reference/RWref_HiC_min1Mb.fasta

samtools dict reference/RWref_HiC_min1Mb.fasta -o reference/RWref_HiC_min1Mb.dict

