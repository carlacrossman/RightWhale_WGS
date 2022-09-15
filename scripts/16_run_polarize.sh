#!/bin/bash
#SBATCH --job-name=polarize_vcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=2:00:00

~/scratch/PY_ENV/polarizeVCF.py --vcf ~/scratch/polarize/all_on_narw_unrelated_filtered_nohaps.vcf.gz --keep outgroup --miss 0.8 -r > all_on_narw_polarized_body.vcf

