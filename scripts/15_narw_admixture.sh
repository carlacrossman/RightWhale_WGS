#!/bin/bash
#SBATCH --job-name=narw_admixture
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=22
#SBATCH --mem=10G
#SBATCH --time=6:00:00

prefix=narw_plink_r.2

cp ${prefix}* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}


for r in {1..10}; do for K in {1..5};
do
	admixture --cv -s ${RANDOM} ${prefix}.bed $K -j20
	mv ${prefix}.${K}.Q ${prefix}.K${K}r${r}.Q
done; done

cp ${prefix}* ${SLURM_SUBMIT_DIR}/out
