#!/bin/bash
# renamebedfiles.sh
# Author: Chris Carson
# Description: Renames ChIP-seq BED files based on experiment target metadata.
# Reads TSV metadata to match file IDs with experiment targets.

# Configuration: Set cell line and paths relative to working directory
CELL_LINE=" "
INPUT_FILE="${CELL_LINE}/metadata.tsv"
FILE_DIR="${CELL_LINE}"

for FILE_PATH in "$FILE_DIR"/*.bed.gz; do
    FILE_NAME=$(basename "$FILE_PATH" .bed.gz)
    
    EXPERIMENT_TARGET=$(awk -v file="$FILE_NAME" -F'\t' '$1 == file {print $23}' "$INPUT_FILE")
    echo "$EXPERIMENT_TARGET"
    
    if [[ -n "$EXPERIMENT_TARGET" ]]; then
        CLEAN_TARGET=${EXPERIMENT_TARGET//-human/}
        NEW_FILE_NAME="$FILE_DIR/${CLEAN_TARGET}-${FILE_NAME}.bed.gz"
        mv "$FILE_PATH" "$NEW_FILE_NAME"
        echo "Renamed $FILE_PATH to $NEW_FILE_NAME"
    else
        echo "No match found for $FILE_NAME in $INPUT_FILE"
    fi
done
