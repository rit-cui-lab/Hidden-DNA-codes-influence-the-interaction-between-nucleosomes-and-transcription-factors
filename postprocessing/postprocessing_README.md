# Postprocessing Scripts Summary

## Overview
Generates GC content profiles around transcription factor binding sites using the outputs of the dnps pipeline. Converts dyad positions to genomic windows, extracts sequences, computes GC content per window, and produces publication-quality plots.

---

## Script List

### Python Scripts

| File | Description |
|------|-------------|
| **calculate_GC_content.py** | Calculates GC content percentage from a FASTA file. Counts G and C nucleotides across all sequences and outputs a single value. |
| **generate_10bp_bed.py** | Converts dyad positions to BED format with 10 bp windows (±5 bp around each dyad). Used for focused GC analysis around TF centers. |
| **generate_147bp_bed.py** | Converts dyad positions to BED format with 147 bp nucleosome windows (±73–74 bp). Duplicate of the preprocessing version; kept here so postprocessing is self-contained. |
| **plotgc.py** | Creates publication-quality GC content plots. Applies 7-element sliding window smoothing and generates symmetrized profiles. |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **generategcscripts.sh** | Generates 200 per-bin SLURM job scripts for each configured protein. Each generated job runs the BED→FASTA→GC pipeline on one bin. |

---

## Pipeline / Data Flow

```
dnps Part 2 output: final_best_site_sorted_unique_{1..200}.bed
    ↓
generate_10bp_bed.py           (convert to 10 bp windows around dyad)
    ↓
bed2fasta                      (extract sequences from hg38)
    ↓
calculate_GC_content.py        (GC% per window)
    ↓
200 per-bin GC files
    ↓
plotgc.py                      (smoothed, symmetrized profile plot)
```

---

## Usage

### Generate per-protein GC analysis jobs

```bash
# 1. Edit CELL_LINE and PROTEINS_TO_PROCESS at the top of generategcscripts.sh
bash generategcscripts.sh

# 2. Submit the generated jobs
for script in ./chipseq/dNPS/<CELL_LINE>/<PROTEIN>_*/job_scripts/*.sh; do
    sbatch "$script"
done

# 3. Plot results for one TF
python plotgc.py <PROTEIN>_GC_aggregated.txt <output.png> <PROTEIN>
```

`generategcscripts.sh` creates 200 SLURM scripts per protein in `<protein_dir>/job_scripts/`. Each script handles one of the 200 bins, so a full protein analysis submits 200 jobs. Aggregate the per-bin GC values into a single file before running `plotgc.py`.

Edit the `#SBATCH` account/partition/email lines in `generategcscripts.sh` for your cluster before submitting.

---

## Inputs and Outputs

**Inputs**
- `final_best_site_sorted_unique_{1..200}.bed` per protein (from `dnps/` Part 2)
- Reference genome FASTA (`hg38.fa`)

**Outputs**
- `<protein>_best_site_sorted_unique_{1..200}_GC.bed.fa` — per-bin FASTA sequences
- `<protein>_best_site_sorted_unique_{1..200}_GC.txt` — per-bin GC fractions
- `<protein>_GC_plot.png` — smoothed symmetric GC profile (via `plotgc.py`)

---

## Dependencies

- Python ≥ 3.9 with numpy, matplotlib
- bed2fasta 5.5.4
- SLURM (for batch submission)
- Spack (optional, for environment management)

See the top-level `README.md` for tested versions and installation.
