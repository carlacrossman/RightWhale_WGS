############################################
###---   Unpacking Bowhead SRA file   ---###
############################################
#
# This script unpacks the bowhead SRA file.
# More information about the file and project can 
# be found on the bowhead genome project website:
# www.bowhead-whale.org
#
# The code below unpacks the SRA file into 6 fastq files:
# Forward and Reverse reads for 3 lanes L003, L004 and L006.
#
# It saves them to have similar file name structures to our data files 
# and then gzips the fastq files.
# 
# Requirements:
#	- SRAToolKit (module: mugqic/sratoolkit/2.8.2-1)
#
##########################################

# Download SRA file

wget 	https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1685383/SRR1685383

# Unpack SRA file
fastq-dump --split-files -v SRR1685383.sra

# Split SRA into fastq files by lane number
awk 'BEGIN {FS = ":"} {lane=$4 ; print > "SRR1685383.L00"lane"_R1_001.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "SRR1685383.L00"lane"_R1_001.fastq"}}' < SRR1685383_1.fastq

awk 'BEGIN {FS = ":"} {lane=$4 ; print > "SRR1685383.L00"lane"_R2_001.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "SRR1685383.L00"lane"_R2_001.fastq"}}' < SRR1685383_2.fastq

#Gzip .fastq files to .gzip files.
gzip SRR1685383.L003_R1_001.fastq
gzip SRR1685383.L003_R2_001.fastq
gzip SRR1685383.L004_R1_001.fastq
gzip SRR1685383.L004_R2_001.fastq
gzip SRR1685383.L006_R1_001.fastq
gzip SRR1685383.L006_R2_001.fastq
