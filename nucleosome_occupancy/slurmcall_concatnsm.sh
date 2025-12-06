#!/bin/bash -l
# slurmcall_concatnsm.sh
# Author: Chris Carson
# Description: Concatenates chromosome-specific occupancy files into single file.
# Combines chr1-22, chrX, and chrY results for BigWig conversion.

#SBATCH -J concat-nsm
#SBATCH -t 7:0:0
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200g
#SBATCH -o concat_output.txt
#SBATCH -e concat_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Configuration: Set cell line to process
CELL_LINE=" "
WORK_DIR="./mnase/${CELL_LINE}"

# Verify working directory exists
if [ ! -d "$WORK_DIR" ]; then
    echo "Error: Working directory $WORK_DIR does not exist"
    exit 1
fi

cd "$WORK_DIR"

OUTPUT_FILE="${CELL_LINE}.bg"
FILES_CONCATENATED=0

echo "Starting concatenation for $CELL_LINE..."
echo "Output file: $OUTPUT_FILE"

# Concatenate autosomal chromosomes (chr1-22)
for i in {1..22}; do
    CHR_FILE="split_filtered_${CELL_LINE}_merged_sorted_dyads_final_chr${i}.bg"
    if [ -f "$CHR_FILE" ]; then
        echo "Adding: $CHR_FILE"
        cat "$CHR_FILE" >> "${OUTPUT_FILE}"
        ((FILES_CONCATENATED++))
    fi
done

# Add chrX if exists
CHR_X_FILE="split_filtered_${CELL_LINE}_merged_sorted_dyads_final_chrX.bg"
if [ -f "$CHR_X_FILE" ]; then
    echo "Adding: $CHR_X_FILE"
    cat "$CHR_X_FILE" >> "${OUTPUT_FILE}"
    ((FILES_CONCATENATED++))
fi

# Add chrY if exists
CHR_Y_FILE="split_filtered_${CELL_LINE}_merged_sorted_dyads_final_chrY.bg"
if [ -f "$CHR_Y_FILE" ]; then
    echo "Adding: $CHR_Y_FILE"
    cat "$CHR_Y_FILE" >> "${OUTPUT_FILE}"
    ((FILES_CONCATENATED++))
fi

echo "Concatenated $FILES_CONCATENATED files"

# Compress intermediate files to save space
echo "Compressing intermediate chromosome files..."
gzip split_filtered_* 2>/dev/null

echo "Concatenation and compression complete!"
echo "Final output: $OUTPUT_FILE"
echo "Compressed files: split_filtered_*.bg.gz"
