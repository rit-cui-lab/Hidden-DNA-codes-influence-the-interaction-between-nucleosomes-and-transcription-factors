# DAC Scripts Summary

## Overview
Scripts for computing the Distance Auto-Correlation (DAC) function of nucleosome dyad positions. DAC quantifies the periodic spacing of nucleosomes along DNA by measuring how often pairs of dyads occur at each distance (0–1000 bp) from each other.

---

## Script List

### 1. **auto_correlation_multiple_new4_corrected_new.pl**
Perl script that computes the distance auto-correlation function. Reads two dyad count files (the same file is passed twice for auto-correlation; two different files for cross-correlation). For every pair of dyads within 1000 bp of each other, multiplies their occurrence counts and accumulates the product in a distance histogram. Outputs a two-column file: distance (bp) and auto-correlation value.

### 2. **my_prog_auto_correlation_1000.sh**
SLURM job submission script. Runs the Perl script sequentially across 24 human chromosomes (chr1–chr22, chrX) using per-chromosome dyad count files. Outputs one auto-correlation file per chromosome.

---

## Data Flow

```
Dyad count files (per chromosome, from preprocessing/)
    ↓
auto_correlation_multiple_new4_corrected_new.pl
    ↓
Per-chromosome auto-correlation files (distance × count)
```

## Usage

### Single file (manual run)

```bash
perl auto_correlation_multiple_new4_corrected_new.pl \
    <input_dyad_file.txt> \
    <input_dyad_file.txt> \
    <output_auto_correlation.txt>
```

Passing the same file twice computes auto-correlation. Passing two different files computes cross-correlation.

**Input format:** tab-separated, three columns — chromosome, dyad position, count (the same format produced by `preprocessing/count_dyads.py`).

**Output format:** tab-separated, two columns — distance (bp, 0–1000), summed correlation count.

### Full chromosome batch (SLURM)

```bash
sbatch my_prog_auto_correlation_1000.sh
```

Edit the `#SBATCH` account, partition, and email lines at the top for your cluster before submitting. The script assumes dyad input files are named `chr{1..22,X}_human_inVitro_sorted_dyad_147bp_type4_new_no_0.txt` in the working directory — adjust the file pattern for your data.

## Key Dependencies

- Perl ≥ 5.26 (tested with 5.26.2)
- SLURM (for batch submission)
- Spack (optional, for environment management)
