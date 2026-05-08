# DAC Scripts Summary

## Overview
Scripts for computing the Distance Auto-Correlation (DAC) function of nucleosome dyad positions. DAC quantifies the periodic spacing of nucleosomes along DNA by measuring how often pairs of dyads occur at each distance (0–1000 bp) from each other.

---

## Script List

### Perl Scripts

| File | Description |
|------|-------------|
| **auto_correlation_multiple_new4_corrected_new.pl** | Computes the distance auto-correlation function. Reads two dyad count files (the same file twice for auto-correlation; two different files for cross-correlation). |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **my_prog_auto_correlation_1000.sh** | SLURM job script. Runs the Perl script across 24 human chromosomes (chr1–22, chrX). |

---

## Pipeline / Data Flow

```
Dyad count files (per chromosome, from preprocessing/)
    ↓
auto_correlation_multiple_new4_corrected_new.pl
    ↓
Per-chromosome auto-correlation files (distance × count)
```

---

## Inputs

### `auto_correlation_multiple_new4_corrected_new.pl`
Three-column TSV: chromosome, dyad position, count (output of `preprocessing/count_dyads.py`). Pass twice for auto-correlation, or pass two different files for cross-correlation. Real example: [`demo/input/HepG2_dyads_chr22.txt`](../demo/input/HepG2_dyads_chr22.txt).

```
chr22	10510092	2
chr22	10510213	1
chr22	10510282	1
```

### `my_prog_auto_correlation_1000.sh`
Per-chromosome dyad count files in the working directory, named `chr{1..22,X}_human_inVitro_sorted_dyad_147bp_type4_new_no_0.txt` (adjust the pattern in the script for your data).

---

## Outputs

### `auto_correlation_multiple_new4_corrected_new.pl`
Two-column TSV (path passed as third argument): distance in bp (0–1000), summed correlation count. Distance 0 dominates (self-pairs); biologically meaningful peaks appear near 147, 294, and 441 bp (nucleosome repeat length and multiples). Real example: [`demo/expected/CTCF_HepG2_chr22_DAC.txt`](../demo/expected/CTCF_HepG2_chr22_DAC.txt).

```
0	28624
1	4260
2	3789
3	3590
4	3453
```

### `my_prog_auto_correlation_1000.sh`
One auto-correlation file per chromosome, named `chr{N}_..._auto_correlation.txt`. Same two-column TSV format as above.

---

## Usage

### Single chromosome (manual)

```bash
perl DAC/auto_correlation_multiple_new4_corrected_new.pl \
    chr22_dyads.txt \
    chr22_dyads.txt \
    chr22_auto_correlation.txt
```

### Full genome (SLURM)

```bash
sbatch DAC/my_prog_auto_correlation_1000.sh
```

Edit the `#SBATCH` account, partition, and email lines before submitting.

### Demo

```bash
bash DAC/run_demo.sh demo/input demo/output/DAC
```

---

## Dependencies

- Perl ≥ 5.26 (tested with 5.26.2)
- SLURM (for batch submission)
- Spack (optional, for environment management)

See the top-level `README.md` for tested versions and installation.
