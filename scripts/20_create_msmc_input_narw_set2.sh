#!/bin/bash
#SBATCH --job-name=generate_MSMC_input_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=21
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

~/scratch/PY_ENV/generate_multihetsep.py \
--mask=EGL254-1_${CHR}_msmcmask.bed.gz \
--mask=EGL272-1_${CHR}_msmcmask.bed.gz \
--mask=EGL140-1_${CHR}_msmcmask.bed.gz \
./EGL254-1_${CHR}.filtered.recode.vcf.gz \
./EGL272-1_${CHR}.filtered.recode.vcf.gz \
./EGL140-1_${CHR}.filtered.recode.vcf.gz > ../input_files/set2_narw_${CHR}_multihetsep.txt
