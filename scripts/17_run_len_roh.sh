#!/bin/bash
#SBATCH --job-name=roh_len
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=24:00:00

python len_roh_python.py
