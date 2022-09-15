#!/bin/bash
#SBATCH --job-name=genotype_gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=24:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_on_narw_gendb_inputs)

FOLDER=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenotypeGVCFs \
   -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
   -V gendb://$FOLDER/ \
   -O ~/scratch/vcf/all_on_narw_HiC_scaffold_$SCAFFOLD.vcf.gz \
   2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/all_on_narw_HiC_scaffold_$SCAFFOLD-gvcf.out
