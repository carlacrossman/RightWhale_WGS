#!/bin/bash
#SBATCH --job-name=IBDSeq
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=3:00:00


##ALSO RAN THIS USING DIFFERENT r2max

vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered.vcf.gz.recode.vcf --max-missing 1.0 --recode --stdout | gzip -c > narw_all_filters_no_missing.vcf.gz

for SCAFFOLD in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=narw_all_filters_no_missing_aug22.vcf.gz out=r1.0/narw_ibdseq_scaff${SCAFFOLD}_1.0 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=1.0;
done

for SCAFFOLD in 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=srw_all_filters_no_missing_aug22.vcf.gz out=r1.0/srw_ibdseq_scaff${SCAFFOLD}_1.0 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=1.0;
done
