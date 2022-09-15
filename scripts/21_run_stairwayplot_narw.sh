#!/bin/bash
#SBATCH --job-name=stariway_plot_unfoldednarw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=12:00:00

bash NARW_unfolded_Aug18.blueprint.sh
