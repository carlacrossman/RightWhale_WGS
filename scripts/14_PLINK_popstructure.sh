#!bin/bash

module load plink/1.9b_6.21-x86_64

# NARW
plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--out narw_plink_r.2

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_plink_r.2.prune.in \
--make-bed --pca --out narw_plink_r.2

# SRW
plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out srw_plink_r.1

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract srw_plink_r.1.prune.in \
--make-bed --pca --out srw_plink_r.1

#ALL

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--remove SRR1685383 \
--out all_plink_r.2

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract all_plink_r.2.prune.in \
--remove SRR1685383 \
--make-bed --pca --out all_plink_r.2
