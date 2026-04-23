#!/bin/bash
# nucleosome_occupancy/run_demo.sh
#
# Demonstrates Gaussian-smoothed nucleosome occupancy calculation on a
# small chr22 slice of HepG2 dyad counts. No Slurm, no cluster required.
#
# Usage:
#   bash nucleosome_occupancy/run_demo.sh <demo_input_dir> <output_dir>
#
# Example:
#   bash nucleosome_occupancy/run_demo.sh demo/input demo/output/nucleosome_occupancy

set -euo pipefail

INPUT_DIR="${1:-demo/input}"
OUTPUT_DIR="${2:-demo/output/nucleosome_occupancy}"

# Resolve all paths to absolute BEFORE any cd, so paths stay valid across cd calls.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTPUT_DIR"
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

# Subsample to keep runtime on a laptop reasonable.
echo "[demo] Subsampling chr22 dyads to first 50k positions..."
head -50000 "$INPUT_DIR/HepG2_dyads_chr22.txt" > "$OUTPUT_DIR/HepG2_dyads_chr22_demo.txt"

echo "[demo] Running calculatensm_multithread.py..."
cd "$OUTPUT_DIR"
python "$SCRIPT_DIR/calculatensm_multithread.py" HepG2_dyads_chr22_demo.txt
cd - > /dev/null

echo "[demo] Normalizing by average..."
python "$SCRIPT_DIR/calculate_occupancy_avg_and_normalize.py" \
    "$OUTPUT_DIR/HepG2_dyads_chr22_demo.bg" \
    "$OUTPUT_DIR/HepG2_dyads_chr22_demo_normalized.bg"

echo ""
echo "[demo] Done. Outputs:"
echo "  $OUTPUT_DIR/HepG2_dyads_chr22_demo.bg             (raw occupancy)"
echo "  $OUTPUT_DIR/HepG2_dyads_chr22_demo_normalized.bg  (normalized)"
echo ""
echo "First 5 lines of normalized output:"
head -5 "$OUTPUT_DIR/HepG2_dyads_chr22_demo_normalized.bg"
