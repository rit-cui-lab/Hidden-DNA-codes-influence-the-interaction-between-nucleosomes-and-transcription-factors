# Postprocessing Scripts Summary

## Overview
Generates GC content profiles around transcription factor binding sites using the outputs of the dnps pipeline. Converts dyad positions to genomic windows, extracts sequences, computes GC content per window, and produces publication-quality plots.

---

## Script List

### Python Scripts

| File | Description |
|------|-------------|
| **calculate_GC_content.py** | Calculates GC content percentage from a FASTA file. |
| **generate_10bp_bed.py** | Converts dyad positions to BED format with 10 bp windows (±5 bp). |
| **generate_147bp_bed.py** | Converts dyad positions to BED format with 147 bp nucleosome windows. Duplicate of the preprocessing version, kept here so the folder is self-contained. |
| **plotgc.py** | Creates publication-quality GC content plots with smoothing and symmetrization. |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **generategcscripts.sh** | Generates 200 per-bin SLURM job scripts for each configured protein. Each script runs the BED→FASTA→GC pipeline on one bin. |

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

## Inputs

### `generate_10bp_bed.py` / `generate_147bp_bed.py`
Three-column BED of dyad positions (output of dnps Part 2). A real demo set is at [`demo/input/dnps_intermediates/`](../demo/input/dnps_intermediates/).

```
chr22	10510092	10510092
chr22	10510213	10510213
```

### `calculate_GC_content.py`
A FASTA file produced by `bed2fasta` from the 10 bp window BED.

```
>chr22:10510087-10510097
ACGTACGTAC
```

### `plotgc.py`
A plain-text file with one GC fraction per line (200 lines for a full ±1000 bp profile).

```
0.42
0.43
0.45
```

### `generategcscripts.sh`
A configured `CELL_LINE` and `PROTEINS_TO_PROCESS` list at the top of the script, plus completed dnps Part 2 output for each protein.

---

## Outputs

### `generate_10bp_bed.py`
A BED file (`<input>_10bp.bed`) where each row spans 10 bp around the dyad.

```
chr22	10510087	10510097
chr22	10510208	10510218
```

### `generate_147bp_bed.py`
A BED file (`<input>_147bp.bed`) where each row spans 147 bp around the dyad.

```
chr22	10510018	10510165
```

### `calculate_GC_content.py`
A single-line file with a single GC fraction (passed as second argument).

```
0.43
```

### `plotgc.py`
A PNG plot of GC content vs distance from TF center, with sliding-window smoothing and symmetrization. Real example: [`demo/expected/CTCF_chr22_GC.png`](../demo/expected/CTCF_chr22_GC.png).

### `generategcscripts.sh`
200 SLURM job scripts in `<protein_dir>/job_scripts/` per configured protein. Each script runs the full BED→FASTA→GC pipeline for one bin.

---

## Usage

### Generate per-protein GC analysis jobs

```bash
# 1. Edit CELL_LINE and PROTEINS_TO_PROCESS at the top of generategcscripts.sh
bash postprocessing/generategcscripts.sh

# 2. Submit the generated jobs
for script in ./chipseq/dNPS/<CELL_LINE>/<PROTEIN>_*/job_scripts/*.sh; do
    sbatch "$script"
done

# 3. Plot results for one TF
python postprocessing/plotgc.py <PROTEIN>_GC_aggregated.txt <output.png> <PROTEIN>
```

### Demo

```bash
bash postprocessing/run_demo.sh demo/input demo/output/postprocessing
```

---

## Dependencies

- Python ≥ 3.9 with numpy, matplotlib
- bed2fasta 5.5.4
- SLURM (for batch submission)
- Spack (optional, for environment management)

See the top-level `README.md` for tested versions and installation.
