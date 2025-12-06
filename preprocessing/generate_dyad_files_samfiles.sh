#!/bin/bash -l
# generate_dyad_files_samfiles.sh
# Author: Chris Carson
# Description: Extracts dyad positions from paired-end SAM files.
# Calculates dyad coordinates from insert size and alignment start position.

#SBATCH -J Dyad-Extraction
#SBATCH -t 7:0:0
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20g
#SBATCH -o output.txt
#SBATCH -e error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Configuration: Set cell line and input file
CELL_LINE=" "
INPUT_FILE="${CELL_LINE}/filtered_${CELL_LINE}.sam"

dir_path=$(dirname "$INPUT_FILE")
filename=$(basename "$INPUT_FILE")
prefix=${filename%.*}
output_file="${dir_path}/${prefix}_dyads.txt"

awk -F'\t' '
/^@/ { next }
{
    value = $9
    if (value > 0) {
        dyad = $4 + int((value - 1) / 2)
        print $3 "\t" dyad
    }
}' "$INPUT_FILE" > "$output_file"
