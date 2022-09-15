#!/bin/bash
#SBATCH --job-name=Combine-MSMC-IM
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-10
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=00:30:00

SET=${SLURM_ARRAY_TASK_ID}

MU="0.9664e-8"
OUT=./rw_msmc_im_set${SET}


~/scratch/PY_ENV/msmc_tools/combineCrossCoal.py ALL_BETWEEN_msmc2_set${SET}.final.txt ALL_NARW_msmc2_set${SET}.final.txt ALL_SRW_msmc2_set${SET}.final.txt > RW_combined_msmc_set${SET}.final.txt
~/scratch/PY_ENV/msmc_tools/MSMC_IM.py -mu $MU -o $OUT --xlog --ylog RW_combined_msmc_set${SET}.final.txt

