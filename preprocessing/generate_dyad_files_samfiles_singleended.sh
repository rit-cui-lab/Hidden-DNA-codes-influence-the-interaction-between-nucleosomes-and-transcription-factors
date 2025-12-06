#!/bin/bash
# generate_dyad_files_samfiles_singleended.sh
# Author: Chris Carson
# Description: Extracts dyad positions from single-ended SAM files.
# Handles forward and reverse strand reads with position adjustments.

#SBATCH -J Dyad-Extraction-SE
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

# Configuration: Set cell line and input/output files
CELL_LINE=" "
INPUT="${CELL_LINE}/${CELL_LINE}.sam"
OUTPUT="${CELL_LINE}/${CELL_LINE}_dyads.txt"

awk -F'\t' -v OFS='\t' '
/^@/ { next }
{
    chr = $3
    flag = $2
    start = $4

    if (flag == 0) {
        dyad = start + 73
    } else if (flag == 16) {
        dyad = start - 39
    } else {
        next
    }

    name = chr "_" dyad
    print chr, dyad-1, dyad, name
}' "$INPUT" > "$OUTPUT"
