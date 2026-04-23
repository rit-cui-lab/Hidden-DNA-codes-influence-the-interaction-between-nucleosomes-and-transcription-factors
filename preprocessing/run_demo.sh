#!/bin/bash
# preprocessing/run_demo.sh
#
# Demonstrates count_dyads.py on a small chr22 slice of HepG2 dyad positions.
# No Slurm, no cluster required.
#
# Usage:
#   bash preprocessing/run_demo.sh <demo_input_dir> <output_dir>
#
# Example:
#   bash preprocessing/run_demo.sh demo/input demo/output/preprocessing

set -euo pipefail

INPUT_DIR="${1:-demo/input}"
OUTPUT_DIR="${2:-demo/output/preprocessing}"

# Resolve all paths to absolute BEFORE any cd.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTPUT_DIR"
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

# The demo dyad file is already counted. To demonstrate count_dyads.py, we
# first expand the counts back into individual dyad occurrences, then re-count.
echo "[demo] Preparing a raw dyad file (one line per occurrence) for count_dyads.py..."
awk 'BEGIN{OFS="\t"} {for (i=0; i<$3; i++) print $1, $2}' \
    "$INPUT_DIR/HepG2_dyads_chr22.txt" \
    > "$OUTPUT_DIR/HepG2_dyads_chr22_raw.txt"

echo "[demo] Running count_dyads.py..."
python "$SCRIPT_DIR/count_dyads.py" "$OUTPUT_DIR/HepG2_dyads_chr22_raw.txt"
# count_dyads.py writes output alongside input, replacing .txt with _final.txt
mv "$OUTPUT_DIR/HepG2_dyads_chr22_raw_final.txt" \
   "$OUTPUT_DIR/HepG2_dyads_chr22_counts.txt"

echo ""
echo "[demo] Done. Output:"
echo "  $OUTPUT_DIR/HepG2_dyads_chr22_counts.txt"
echo ""
echo "First 5 lines:"
head -5 "$OUTPUT_DIR/HepG2_dyads_chr22_counts.txt"
echo ""
echo "Expected to match demo/expected/HepG2_dyads_chr22_counts.head.txt (first 100 lines)."
