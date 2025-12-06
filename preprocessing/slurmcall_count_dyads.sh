#!/bin/bash -l
# slurmcall_count_dyads.sh
# Author: Chris Carson
# Description: SLURM job submission wrapper for count_dyads.py.
# Configures computational resources and calls the Python script.

#SBATCH -J Dyad-Count
#SBATCH -t 1-00:00:00
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=350g
#SBATCH -o output.txt
#SBATCH -e error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

CELL_LINE=" "
INPUT_FILE="${CELL_LINE}/filtered_${CELL_LINE}_merged_sorted_dyads.txt"

cd "./${CELL_LINE}/"
python ./scripts/preprocessing/count_dyads.py "$INPUT_FILE"
