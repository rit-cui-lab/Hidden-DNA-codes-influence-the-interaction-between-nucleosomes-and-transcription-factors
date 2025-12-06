#!/bin/bash
# splitdyadfiles.sh
# Author: Chris Carson
# Description: Splits large dyad position files by chromosome.
# Creates separate files for chr1-22, chrX, and chrY for parallel processing.

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    echo "Example: $0 ./mnase/HepG2/filtered_HepG2_dyads_final.txt"
    exit 1
fi

INPUT_FILE="$1"
INPUT_DIR=$(dirname "$INPUT_FILE")
BASE_FILENAME=$(basename "$INPUT_FILE" .txt)
PREFIX="${INPUT_DIR}/split_${BASE_FILENAME}"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist"
    exit 1
fi

echo "Splitting $INPUT_FILE by chromosome..."
echo "Output prefix: $PREFIX"

# Array of chromosomes to process
CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" \
             "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" \
             "chr20" "chr21" "chr22" "chrX" "chrY")

FILE_COUNT=0

# Process each chromosome
for CHR in "${CHROMOSOMES[@]}"; do
    OUTPUT_FILE="${PREFIX}_${CHR}.txt"
    
    # Filter lines matching this chromosome
    awk -v chr="$CHR" '$1 == chr {print}' "$INPUT_FILE" > "$OUTPUT_FILE"
    
    # Remove empty files and count created files
    if [ ! -s "$OUTPUT_FILE" ]; then
        rm "$OUTPUT_FILE"
    else
        echo "Created: $OUTPUT_FILE"
        ((FILE_COUNT++))
    fi
done

echo "Processing complete!"
echo "Created $FILE_COUNT chromosome files"
