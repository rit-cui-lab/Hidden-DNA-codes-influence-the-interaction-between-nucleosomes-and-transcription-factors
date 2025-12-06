#!/bin/bash -l
# filtersamfiles.sh
# Author: Chris Carson
# Description: Filters SAM files based on insert size range (120-180 bp).
# Converts to BAM format, sorts, and indexes. Removes original SAM files.

#SBATCH -J SAM-Filter
#SBATCH -t 2-7:0:0
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=75g
#SBATCH -o output.txt
#SBATCH -e error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load samtools

# Configuration: Set cell line and paths relative to working directory
CELL_LINE="HepG2"
INPUT_FILE="${CELL_LINE}/${CELL_LINE}.sam"

input_dir=$(dirname "$INPUT_FILE")
input_name=$(basename "$INPUT_FILE" .sam)
filtered_sam="${input_dir}/filtered_${input_name}.sam"

echo "Converting original SAM to BAM..."
original_bam="${input_dir}/${input_name}.bam"
samtools view -bS "$INPUT_FILE" > "$original_bam"
echo "Original BAM file created: $original_bam"

echo "Filtering SAM file..."
awk 'BEGIN {OFS="\t"}
    /^@/ {print; next}
    $9 >= 120 && $9 <= 180 || $9 <= -120 && $9 >= -180 {print}
    ' "$INPUT_FILE" > "$filtered_sam"
echo "Filtered SAM file created: $filtered_sam"

echo "Sorting original BAM file..."
samtools sort "$original_bam" -o "${input_dir}/${input_name}_sorted.bam"

echo "Indexing sorted original BAM file..."
samtools index "${input_dir}/${input_name}_sorted.bam"

echo "Removing original SAM file to save space..."
rm "$INPUT_FILE"

echo "Process completed!"
echo "Files created:"
echo "- Original BAM: $original_bam"
echo "- Sorted original BAM: ${input_dir}/${input_name}_sorted.bam"
echo "- BAM index: ${input_dir}/${input_name}_sorted.bam.bai"
echo "- Filtered SAM: $filtered_sam"
