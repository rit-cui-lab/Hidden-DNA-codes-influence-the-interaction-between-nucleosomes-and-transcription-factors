#!/bin/bash
# postprocessing/run_demo.sh
#
# Demonstrates GC content computation and plotting across bins 90-110 around
# CTCF binding site centers (chr22, HepG2). Runs end-to-end without Slurm.
#
# PREREQUISITES (must be on PATH before running this script):
#   - python with numpy and matplotlib
#   - bed2fasta (tested with 5.5.4)
#
# On the RIT SPORC cluster this is provided by:
#   spack env activate factor-x86_64-24062401
#
# On other systems, install bed2fasta separately and ensure it is on PATH.
#
# Usage:
#   bash postprocessing/run_demo.sh <demo_input_dir> <output_dir>
#
# Example:
#   bash postprocessing/run_demo.sh demo/input demo/output/postprocessing

set -euo pipefail

# Prerequisite check: fail early with a helpful message if bed2fasta is missing.
if ! command -v bed2fasta > /dev/null 2>&1; then
    echo "ERROR: bed2fasta not found on PATH." >&2
    echo "Install bed2fasta (tested with 5.5.4) or, on the RIT cluster, run:" >&2
    echo "    spack env activate factor-x86_64-24062401" >&2
    exit 1
fi

INPUT_DIR="${1:-demo/input}"
OUTPUT_DIR="${2:-demo/output/postprocessing}"

# Resolve all paths to absolute BEFORE any cd.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTPUT_DIR"
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
GENOME_FA="$INPUT_DIR/chr22.fa"

# For each bin: convert dyads to 10 bp windows, extract FASTA, compute GC.
echo "[demo] Computing GC content for bins 90-110..."
> "$OUTPUT_DIR/CTCF_chr22_GC_values.txt"
for i in $(seq 90 110); do
    bin_bed="$INPUT_DIR/dnps_intermediates/final_best_site_sorted_unique_${i}.bed"

    python "$SCRIPT_DIR/generate_10bp_bed.py" "$bin_bed"
    mv "${bin_bed%.bed}_10bp.bed" "$OUTPUT_DIR/"

    bed2fasta -s -both \
        -o "$OUTPUT_DIR/bin_${i}_GC.bed.fa" \
        "$OUTPUT_DIR/final_best_site_sorted_unique_${i}_10bp.bed" \
        "$GENOME_FA"

    python "$SCRIPT_DIR/calculate_GC_content.py" \
        "$OUTPUT_DIR/bin_${i}_GC.bed.fa" \
        "$OUTPUT_DIR/bin_${i}_GC.txt"

    cat "$OUTPUT_DIR/bin_${i}_GC.txt" >> "$OUTPUT_DIR/CTCF_chr22_GC_values.txt"
done

# plotgc.py expects 200 values; pad with edge values so the smoothing window
# has something to chew on.
echo "[demo] Padding 21 bin values to 200 for plotgc.py..."
python - "$OUTPUT_DIR" <<'PYEOF'
import sys, os
out_dir = sys.argv[1]
with open(os.path.join(out_dir, "CTCF_chr22_GC_values.txt")) as f:
    vals = [float(line.strip()) for line in f if line.strip()]
# Center the 21 real values at bins 90-110 of a 200-element array,
# pad the rest with the edge values (reasonable for a demo).
padded = [vals[0]] * 89 + vals + [vals[-1]] * (200 - 89 - len(vals))
assert len(padded) == 200, len(padded)
with open(os.path.join(out_dir, "CTCF_chr22_GC_padded.txt"), "w") as f:
    for v in padded:
        f.write(f"{v}\n")
PYEOF

echo "[demo] Generating GC plot..."
python "$SCRIPT_DIR/plotgc.py" \
    "$OUTPUT_DIR/CTCF_chr22_GC_padded.txt" \
    "$OUTPUT_DIR/CTCF_chr22_GC.png" \
    "CTCF (demo)"

echo ""
echo "[demo] Done. Output plot:"
echo "  $OUTPUT_DIR/CTCF_chr22_GC.png"
echo ""
echo "Note: the demo covers only bins 90-110 (±100 bp around CTCF center)."
echo "Outside that range the plot shows padded values, not real data."
echo "Compare the central region against demo/expected/CTCF_chr22_GC.png"
