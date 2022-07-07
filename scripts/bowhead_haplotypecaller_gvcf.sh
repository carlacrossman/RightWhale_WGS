#!/bin/bash
#SBATCH --job-name=bowhead-haplotype-gvcf_srw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21%10
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=24:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bowhead_sample_scaffolds_SRW)

SAMPLE=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked_SRW.bam \
  -O ~/scratch/gvcf/bowhead/map_to_SRW/$SCAFFOLD/$SAMPLE-$SCAFFOLD-SRW.g.vcf.gz \
  -ERC GVCF \
  -L $SCAFFOLD \
  2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/$SAMPLE-$SCAFFOLD-hapcall-gvcf_SRW.out
