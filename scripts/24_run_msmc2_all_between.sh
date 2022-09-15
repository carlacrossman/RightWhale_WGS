#!/bin/bash
#SBATCH --job-name=MSMC2_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-10
#SBATCH --cpus-per-task=20
#SBATCH --mem=124G
#SBATCH --time=12:00:00

SET=${SLURM_ARRAY_TASK_ID}

cp input_files/all*set${SET}* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}

msmc2 -t 20 -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -o ALL_BETWEEN_msmc2_set${SET} \
        all_1_multihetsep_set${SET}.txt \
        all_2_multihetsep_set${SET}.txt \
        all_3_multihetsep_set${SET}.txt \
        all_4_multihetsep_set${SET}.txt \
        all_5_multihetsep_set${SET}.txt \
        all_6_multihetsep_set${SET}.txt \
        all_7_multihetsep_set${SET}.txt \
	all_8_multihetsep_set${SET}.txt \
        all_9_multihetsep_set${SET}.txt \
        all_10_multihetsep_set${SET}.txt \
	all_11_multihetsep_set${SET}.txt \
        all_12_multihetsep_set${SET}.txt \
        all_13_multihetsep_set${SET}.txt \
        all_14_multihetsep_set${SET}.txt \
        all_15_multihetsep_set${SET}.txt \
        all_16_multihetsep_set${SET}.txt \
        all_17_multihetsep_set${SET}.txt \
        all_18_multihetsep_set${SET}.txt \
        all_19_multihetsep_set${SET}.txt \
        all_20_multihetsep_set${SET}.txt \
        all_22_multihetsep_set${SET}.txt

cp ALL_BETWEEN_msmc2_set${SET}* ${SLURM_SUBMIT_DIR}/results/
