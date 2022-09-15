#!/bin/bash
#SBATCH --job-name=IBDNE_NARW
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=01:00:00

#cat ~/scratch/IBDNe/IBDSeq/r1.0/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_0.8_aug29.map out=narw_r2max1.0_recom0.8_cm2 nthreads=20 mincm=2

#cat ~/scratch/IBDNe/IBDSeq/r1.0/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_0.9_aug29.map out=narw_r2max1.0_recom0.9_cm2 nthreads=20 mincm=2

#cat ~/scratch/IBDNe/IBDSeq/r1.0/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_1.0_aug29.map out=narw_r2max1.0_recom1.0_cm2 nthreads=20 mincm=2

cat ~/scratch/IBDNe/IBDSeq/r1.0/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_1.1_aug29.map out=narw_r2max1.0_recom1.1_cm2 nthreads=20 mincm=2

cat ~/scratch/IBDNe/IBDSeq/r1.0/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_1.2_aug29.map out=narw_r2max1.0_recom1.2_cm2 nthreads=20 mincm=2

