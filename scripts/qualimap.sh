##########################################
###--- Running Qualimap on all .bams---###
##########################################
#
# This script runs qualimap on all merged .bam files.
# 
# This code uses a file called qualimap_sample_file 
# saved in the merged_bam/
#
# Requirements:
#	- Qualimap (module: mugqic/qualimap/2.2.1)
#	- file 'qualimap_sample_file' that specifies samples, file name and population
#
# To be run from: merged_bam/
#
##########################################

qualimap multi-bamqc -d qualimap_sample_file \
    -outdir QC/mapQC -nr 10000 \
    --java-mem-size=25G -r
