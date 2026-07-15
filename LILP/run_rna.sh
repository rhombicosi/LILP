#!/bin/bash
#SBATCH --job-name=rna_ilp
#SBATCH --partition=cpu_zen2
#SBATCH --array=0-3
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=7-00:00:00

python store_results.py $SLURM_ARRAY_TASK_ID