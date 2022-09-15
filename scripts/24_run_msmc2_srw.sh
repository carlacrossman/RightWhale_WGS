#!/bin/bash
#SBATCH --job-name=MSMC2_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-10
#SBATCH --cpus-per-task=20
#SBATCH --mem=124G
#SBATCH --time=24:00:00

SET=${SLURM_ARRAY_TASK_ID}

cp input_files/all*_set${SET}* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}

msmc2 -t 20 -I 0,1,2,3 -o ALL_NARW_msmc2_set${SET} \
        all_1_multihetsep_set${SET}.txt \
        all_2_multihetsep${SET}.txt \
        all_3_multihetsep${SET}.txt \
        all_4_multihetsep${SET}.txt \
        all_5_multihetsep${SET}.txt \
        all_6_multihetsep${SET}.txt \
        all_7_multihetsep${SET}.txt \
	all_8_multihetsep${SET}.txt \
        all_9_multihetsep${SET}.txt \
        all_10_multihetsep${SET}.txt \
	all_11_multihetsep${SET}.txt \
        all_12_multihetsep${SET}.txt \
        all_13_multihetsep${SET}.txt \
        all_14_multihetsep${SET}.txt \
        all_15_multihetsep${SET}.txt \
        all_16_multihetsep${SET}.txt \
        all_17_multihetsep${SET}.txt \
        all_18_multihetsep${SET}.txt \
        all_19_multihetsep${SET}.txt \
        all_20_multihetsep${SET}.txt \
        all_22_multihetsep${SET}.txt

cp ALL_NARW_msmc2_set${SET} ${SLURM_SUBMIT_DIR}/results/
