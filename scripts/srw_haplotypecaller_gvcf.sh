###########################################
###---Calling haplotypes in GVCF mode---###
###########################################
#
# This script uses GATK4 HaplotypeCaller
# to call haplotypes on all sites in gvcf mode
# in each scaffold independently.  
# Standard Error is saved to QC/gvcfQC/.
# 
# Requirements:
#	- java (module: mugqic/java/openjdk-jdk1.8.0_72)
#	- GATK4 (module: mugqic/GenomeAnalysisTK/4.1.0.0)
#	- file 'srw_sample_scaffolds' that specifies samples by lane to be sorted
#	- merged .bam files in merged_bam/
#	- reference genome saved as reference/RWref_HiC_min1Mb.fasta
#	- This was run as an array job with a job for each scaffold for each individual
#
# To be run from: main project directory
#
##########################################

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_sample_scaffolds)

SAMPLE=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
  -R reference/RWref_HiC_min1Mb.fasta \
  -I merged_bam/$SAMPLE-merged_marked.bam \
  -O ~/scratch/gvcf/SRW/$SCAFFOLD/$SAMPLE-$SCAFFOLD.g.vcf.gz \
  -ERC GVCF \
  -L $SCAFFOLD \
  2> QC/gvcfQC/$SAMPLE-$SCAFFOLD-hapcall-gvcf.out
