#!/bin/bash
#SBATCH --job-name=rw-haplotype-gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=245
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=12:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rw_sample_scaffolds)

SAMPLE=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
  -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam \
  -O ~/scratch/gvcf/NARW/$SCAFFOLD/$SAMPLE-$SCAFFOLD.g.vcf.gz \
  -ERC GVCF \
  -L $SCAFFOLD \
  2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/$SAMPLE-$SCAFFOLD-hapcall-gvcf.out
