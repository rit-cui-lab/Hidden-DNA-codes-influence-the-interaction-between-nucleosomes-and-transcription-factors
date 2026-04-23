#!/bin/bash
# dnps/run_demo.sh
#
# Demonstrates dnps Part 3 (nucleosome 4-type classification) and ΔNPS
# computation on 21 pre-computed bins (90-110, spanning -100 bp to +100 bp
# around CTCF binding site centers on chr22 in HepG2).
#
# Primary output is a tab-separated table of type counts and ΔNPS values
# that reviewers can diff exactly against demo/expected/CTCF_dnps_demo.tsv.
# A plot is also generated for visual inspection.
#
# Runs end-to-end without Slurm.
#
# PREREQUISITES (must be on PATH before running this script):
#   - python with numpy and matplotlib
#   - perl
#   - bed2fasta (tested with 5.5.4)
#
# On the RIT SPORC cluster this is provided by:
#   spack env activate factor-x86_64-24062401
#
# On other systems, install bed2fasta separately and ensure it is on PATH.
#
# Usage:
#   bash dnps/run_demo.sh <demo_input_dir> <output_dir>
#
# Example:
#   bash dnps/run_demo.sh demo/input demo/output/dnps

set -euo pipefail

# Prerequisite check: fail early with a helpful message if bed2fasta is missing.
if ! command -v bed2fasta > /dev/null 2>&1; then
    echo "ERROR: bed2fasta not found on PATH." >&2
    echo "Install bed2fasta (tested with 5.5.4) or, on the RIT cluster, run:" >&2
    echo "    spack env activate factor-x86_64-24062401" >&2
    exit 1
fi

INPUT_DIR="${1:-demo/input}"
OUTPUT_DIR="${2:-demo/output/dnps}"

# Resolve all paths to absolute BEFORE any cd.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTPUT_DIR"
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
GENOME_FA="$INPUT_DIR/chr22.fa"

# Part 3: for each bin, expand dyads to 147 bp footprints, extract sequences,
# classify into 4 types based on WW/SS dinucleotide patterns.
echo "[demo] Running Part 3 (classify nucleosomes into 4 types) for bins 90-110..."
for i in $(seq 90 110); do
    bin_bed="$INPUT_DIR/dnps_intermediates/final_best_site_sorted_unique_${i}.bed"

    python "$SCRIPT_DIR/generate_147bp_bed.py" "$bin_bed"
    mv "${bin_bed%.bed}_147bp.bed" "$OUTPUT_DIR/"

    bed2fasta -s -both \
        -o "$OUTPUT_DIR/best_site_sorted_unique_${i}_147bp.bed.fa" \
        "$OUTPUT_DIR/final_best_site_sorted_unique_${i}_147bp.bed" \
        "$GENOME_FA"

    perl "$SCRIPT_DIR/separate_nucleosomes.pl" \
        "$OUTPUT_DIR/best_site_sorted_unique_${i}_147bp.bed.fa" \
        "$OUTPUT_DIR/best_site_sorted_unique_${i}_147bp_four_types.txt"
done

# Produce a diffable table + a reference plot.
echo "[demo] Aggregating four-types counts into CTCF_dnps_demo.tsv..."
python - "$OUTPUT_DIR" <<'PYEOF'
import os, sys
import matplotlib.pyplot as plt

out_dir = sys.argv[1]
bins = list(range(90, 111))

rows = []
for i in bins:
    path = os.path.join(out_dir, f"best_site_sorted_unique_{i}_147bp_four_types.txt")
    with open(path) as fh:
        parts = [int(x) for x in fh.readline().strip().split("\t")[:4]]
    total = sum(parts)
    # Bin i (1-indexed in the 200-bin pipeline) corresponds to
    # (-1000 + 10*(i-1)) bp from the TF center.
    distance = -1000 + 10 * (i - 1)
    dnps = (parts[0] - parts[3]) / total * 100 if total > 0 else 0.0
    rows.append((i, distance, parts[0], parts[1], parts[2], parts[3], total, dnps))

# Write a tab-separated, deterministic table. Use fixed precision on ΔNPS
# so reviewers can compare outputs with a plain diff.
table_path = os.path.join(out_dir, "CTCF_dnps_demo.tsv")
with open(table_path, "w") as fh:
    fh.write("bin\tdistance_bp\ttype1\ttype2\ttype3\ttype4\ttotal\tdelta_NPS_pct\n")
    for row in rows:
        fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6f}\n".format(*row))
print(f"[demo] Wrote {table_path}")

# Reference plot (useful for visual inspection, not for diffing).
xs = [r[1] for r in rows]
ys = [r[7] for r in rows]
plt.figure(figsize=(8, 4))
plt.plot(xs, ys, "o-", color="blue")
plt.axhline(0, color="grey", linestyle=":", linewidth=0.8)
plt.axvline(0, color="grey", linestyle=":", linewidth=0.8)
plt.xlabel("Distance from CTCF center (bp)")
plt.ylabel(r"$\Delta$NPS, %")
plt.title("CTCF / HepG2 / chr22 — demo (bins 90–110, unsmoothed)")
plt.tight_layout()
plot_path = os.path.join(out_dir, "CTCF_dnps_demo.png")
plt.savefig(plot_path, dpi=150)
print(f"[demo] Wrote {plot_path}")
PYEOF

echo ""
echo "[demo] Done."
echo ""
echo "Primary output (diffable table):"
echo "  $OUTPUT_DIR/CTCF_dnps_demo.tsv"
echo ""
echo "First 5 rows:"
head -5 "$OUTPUT_DIR/CTCF_dnps_demo.tsv"
echo ""
echo "Reference plot (visual only):"
echo "  $OUTPUT_DIR/CTCF_dnps_demo.png"
echo ""
echo "To verify reproducibility:"
echo "  diff $OUTPUT_DIR/CTCF_dnps_demo.tsv demo/expected/CTCF_dnps_demo.tsv"
echo ""
echo "The plot will be noisy at this sample size (~300 nucleosomes/bin);"
echo "the real pipeline averages across many TFs and thousands of sites."
