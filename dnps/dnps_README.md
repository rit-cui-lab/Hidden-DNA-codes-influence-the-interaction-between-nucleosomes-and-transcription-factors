# dNPS Pipeline Scripts Summary

## Overview
Implements the delta-NPS (ΔNPS, Nucleosome Positioning Score) pipeline for analyzing nucleosome positioning patterns around transcription factor binding sites. For each TF, the pipeline bins the region ±1000 bp around TF centers into 200 × 10 bp windows, identifies the nucleosomes whose dyads fall in each window, classifies those nucleosomes into four dinucleotide pattern types, and computes a directional score (ΔNPS) measuring rotational orientation of nucleosomes relative to the TF.

---

## Script List

### Python Scripts

| File | Description |
|------|-------------|
| **binning_bed.py** | Creates 200 binned BED files from TF binding sites (10 bp windows spanning ±1000 bp). |
| **generate_147bp_bed.py** | Converts dyad positions to 147 bp nucleosome footprints. |
| **makeavgdnpsprofile.py** | Generates averaged profile plots for groups of TFs. Three-axis plot: in vivo (red), ΔNPS (blue), in vitro (green). |
| **makeindividualplots.py** | Generates per-TF profile plots. |

### Perl Scripts

| File | Description |
|------|-------------|
| **separate_nucleosomes.pl** | Classifies each 147 bp nucleosome sequence into one of four types based on WW/SS dinucleotides at hardcoded rotational positions. |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **callpipeline1.sh** | Part 1 — submits SLURM jobs to run `binning_bed.py` for each TF. |
| **callpipeline2.sh** | Part 2 — multi-node SLURM array jobs that intersect binned TF regions with MNase dyads via `bedtools intersect`. |
| **callpipeline3.sh** | Part 3 — converts dyads to 147 bp footprints, extracts sequences with `bed2fasta`, classifies with `separate_nucleosomes.pl`. |
| **makeavgplotscripts.sh** | Submits per-(cell line, profile) SLURM jobs for `makeavgdnpsprofile.py`. |
| **callindividualplots.sh** | Submits per-TF SLURM jobs for `makeindividualplots.py`. |
| **call3lineplots.sh** | Alternative averaged-plot submission. |

---

## Pipeline / Data Flow

```
Per-TF MEME files + in vivo/in vitro BigWigs (from nucleosome_occupancy/)
    ↓
Part 1 — callpipeline1.sh → binning_bed.py        (200 × 10 bp windows per TF)
    ↓
Part 2 — callpipeline2.sh → bedtools intersect    (intersect bins with MNase dyads)
    ↓
Part 3 — callpipeline3.sh →
    ↓   generate_147bp_bed.py   (147 bp footprints)
    ↓   bed2fasta               (extract sequences)
    ↓   separate_nucleosomes.pl (classify into 4 types)
    ↓
200 four-type classification files per TF
    ↓
Visualization:
    ├── makeavgplotscripts.sh   → makeavgdnpsprofile.py  (averaged plots)
    ├── callindividualplots.sh  → makeindividualplots.py (per-TF plots)
    └── call3lineplots.sh       → make3lineplots.py      (alternative plots)
```

---

## Inputs

### `binning_bed.py` (Part 1)
Three-column BED of TF binding sites (`best_site_sorted_unique.bed` from `nucleosome_occupancy/`) plus a TF name and cell line name as command-line arguments.

```
chr22	41403960	41404006
chr22	42091010	42091056
```

### `callpipeline2.sh` (Part 2)
The 200 binned BED files from Part 1, plus per-chromosome MNase dyad bedGraphs (`split_filtered_<cell_line>_dyads_final_chr<N>.bg` from `nucleosome_occupancy/`).

### `generate_147bp_bed.py` (Part 3)
Three-column BED of dyad positions (output of Part 2). A real demo set (bins 90–110) is at [`demo/input/dnps_intermediates/`](../demo/input/dnps_intermediates/).

```
chr22	10510092	10510092
chr22	10510213	10510213
```

### `separate_nucleosomes.pl`
A FASTA file produced by `bed2fasta` from the 147 bp footprint BED. Each record is a 147 bp nucleosome sequence.

```
>chr22:10510018-10510165
ACGTACGTACGT...
```

### `makeavgdnpsprofile.py` / `makeindividualplots.py`
A TF list file (one TF and meme name per line, space-separated) and a TF-to-profile mapping file (`avgs/tf_lists/tf_profile_mapping_<cell_line>.txt`). Plus the 200 four_types files per TF produced by Part 3.

---

## Profile Categories

TFs are classified into five profile types based on the shape of their in vivo vs in vitro occupancy:

| Category    | Meaning                                 |
|-------------|-----------------------------------------|
| **DIP_DIP**  | Dip in both in vivo and in vitro |
| **DIP_PEAK** | Dip in vivo, peak in vitro |
| **PEAK_DIP** | Peak in vivo, dip in vitro |
| **PEAK**     | Peak in both |
| **AMBIGUOUS** | Unclear or mixed |

---

## Usage

```bash
bash dnps/callpipeline1.sh           # binning
bash dnps/callpipeline2.sh           # intersection
bash dnps/callpipeline3.sh           # classification
bash dnps/makeavgplotscripts.sh      # averaged plots
bash dnps/callindividualplots.sh     # per-TF plots
```

Edit the `CONFIGURATION` block at the top of each script (BASE_DIR, EMAIL, ACCOUNT, PARTITION, SPACK_ENV) and `BASE_DIR` inside `makeavgdnpsprofile.py` and `makeindividualplots.py`.

### Demo

```bash
bash dnps/run_demo.sh demo/input demo/output/dnps
```

Demo runs Part 3 only on 21 pre-computed bins (90–110) and produces a diffable TSV (`CTCF_dnps_demo.tsv`) plus a reference plot.

---

## Dependencies

- Python ≥ 3.9 with numpy, matplotlib
- Perl ≥ 5.26
- bedtools ≥ 2.31
- bed2fasta 5.5.4
- SLURM (multi-node array jobs required for Part 2)
- Spack (optional, for environment management)

See the top-level `README.md` for tested versions and installation.
