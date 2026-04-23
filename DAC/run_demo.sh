#!/bin/bash
# DAC/run_demo.sh
#
# Demonstrates Distance Auto-Correlation (DAC) on a chr22 slice of HepG2
# MNase-seq dyad counts. No Slurm, no cluster required.
#
# Usage:
#   bash DAC/run_demo.sh <demo_input_dir> <output_dir>
#
# Example:
#   bash DAC/run_demo.sh demo/input demo/output/DAC

set -euo pipefail

INPUT_DIR="${1:-demo/input}"
OUTPUT_DIR="${2:-demo/output/DAC}"

# Resolve all paths to absolute BEFORE any cd.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTPUT_DIR"
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

# Subsample to keep runtime on a laptop reasonable. The Perl script is O(N^2)
# over pairs within 1000 bp, so runtime scales quickly with input size.
# 20k dyads cover ~2 Mb of chr22 and run in ~1 minute on a laptop.
echo "[demo] Subsampling chr22 dyads to first 20k positions..."
head -20000 "$INPUT_DIR/HepG2_dyads_chr22.txt" \
    > "$OUTPUT_DIR/HepG2_dyads_chr22_demo.txt"

echo "[demo] Running auto-correlation..."
perl "$SCRIPT_DIR/auto_correlation_multiple_new4_corrected_new.pl" \
    "$OUTPUT_DIR/HepG2_dyads_chr22_demo.txt" \
    "$OUTPUT_DIR/HepG2_dyads_chr22_demo.txt" \
    "$OUTPUT_DIR/HepG2_chr22_auto_correlation.txt" > /dev/null

echo ""
echo "[demo] Done. Output:"
echo "  $OUTPUT_DIR/HepG2_chr22_auto_correlation.txt"
echo ""
echo "First 15 lines (distances 0-14 bp):"
head -15 "$OUTPUT_DIR/HepG2_chr22_auto_correlation.txt"
echo ""
echo "Expected: strong peaks near 147 bp and its multiples (nucleosome repeat length)."
