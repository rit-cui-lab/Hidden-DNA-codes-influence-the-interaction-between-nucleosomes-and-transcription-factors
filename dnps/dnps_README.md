# dNPS Pipeline Scripts Summary

## Overview
Implements the delta-NPS (ΔNPS, Nucleosome Positioning Score) pipeline for analyzing nucleosome positioning patterns around transcription factor binding sites. For each TF, the pipeline bins the region ±1000 bp around TF centers into 200 × 10 bp windows, identifies the nucleosomes whose dyads fall in each window, classifies those nucleosomes into four dinucleotide pattern types, and computes a directional score (ΔNPS) that measures the rotational orientation of nucleosomes relative to the TF. Output plots combine in vivo occupancy, in vitro occupancy, and ΔNPS on a single axis.

---

## Script List

### Python Scripts

| File | Description |
|------|-------------|
| **binning_bed.py** | Creates 200 binned BED files from TF binding sites. Each bin is a 10 bp window at a specific offset spanning −1000 bp to +1000 bp around the TF center. Foundation of the spatial analysis. |
| **generate_147bp_bed.py** | Converts dyad positions to 147 bp nucleosome footprints (±73–74 bp around each dyad). |
| **makeavgdnpsprofile.py** | Generates averaged profile plots for groups of TFs in the same profile category. Three-axis plot: in vivo occupancy (red), ΔNPS (blue), in vitro occupancy (green). |
| **makeindividualplots.py** | Generates per-TF profile plots. Handles missing in vitro data gracefully and logs TFs with missing data to a tracking file. |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **callpipeline1.sh** | Part 1 — submits SLURM jobs to run `binning_bed.py` for each TF, producing 200 binned BED files. |
| **callpipeline2.sh** | Part 2 — submits multi-node SLURM array jobs that intersect binned TF regions with MNase-seq dyad positions via `bedtools intersect`. Produces `final_best_site_sorted_unique_{1..200}.bed`. |
| **callpipeline3.sh** | Part 3 — submits SLURM jobs that convert dyads to 147 bp footprints, extract sequences with `bed2fasta`, and classify each nucleosome into one of four dinucleotide pattern types via `separate_nucleosomes.pl`. |
| **makeavgplotscripts.sh** | Submits SLURM jobs that run `makeavgdnpsprofile.py` per (cell line, profile category) combination. |
| **callindividualplots.sh** | Submits SLURM jobs that run `makeindividualplots.py` per TF. |
| **call3lineplots.sh** | Alternative averaged-plot submission that calls `make3lineplots.py` (alternative visualization). |

---

## Pipeline / Data Flow

```
Per-TF MEME files + in vivo/in vitro BigWigs (from nucleosome_occupancy/)
    ↓
Part 1 — callpipeline1.sh → binning_bed.py
    ↓                        (200 × 10 bp windows per TF)
Part 2 — callpipeline2.sh → bedtools intersect
    ↓                        (intersect bins with MNase dyad bedGraphs)
Part 3 — callpipeline3.sh →
    ↓   generate_147bp_bed.py   (147 bp footprints)
    ↓   bed2fasta               (extract sequences)
    ↓   separate_nucleosomes.pl (classify into 4 types)
    ↓
200 four-type classification files per TF
    ↓
Visualization:
    ├── makeavgplotscripts.sh   → makeavgdnpsprofile.py  (grouped avg plots)
    ├── callindividualplots.sh  → makeindividualplots.py (per-TF plots)
    └── call3lineplots.sh       → make3lineplots.py      (alternative plots)
```

---

## Usage

### Typical run order

```bash
# Parts 1–3 run sequentially; each generates and submits SLURM jobs
bash callpipeline1.sh           # binning
bash callpipeline2.sh           # intersection
bash callpipeline3.sh           # classification

# Visualization — run after Parts 1–3 finish for all TFs
bash makeavgplotscripts.sh      # averaged plots by profile category
bash callindividualplots.sh     # per-TF plots
```

Before submitting, edit the `CONFIGURATION` block at the top of each script:

```bash
BASE_DIR="/path/to/your/project"
EMAIL="your_email@example.com"
ACCOUNT="your_account"
PARTITION="your_partition"
SPACK_ENV="your_spack_env"
```

Also edit `BASE_DIR` inside `makeavgdnpsprofile.py` and `makeindividualplots.py`.

---

## Profile Categories

TFs are classified into five profile types based on the shape of their in vivo vs in vitro occupancy around the binding site:

| Category    | Meaning                                 |
|-------------|-----------------------------------------|
| **DIP_DIP**  | Dip in both in vivo and in vitro occupancy |
| **DIP_PEAK** | Dip in vivo, peak in vitro               |
| **PEAK_DIP** | Peak in vivo, dip in vitro               |
| **PEAK**     | Peak in both                             |
| **AMBIGUOUS** | Unclear or mixed pattern                |

Mappings are stored in `avgs/tf_lists/tf_profile_mapping_<cell_line>.txt` (one TF and category per line, tab-separated).

---

## Inputs and Outputs

**Inputs**
- Per-TF MEME files (from `preprocessing/`)
- Per-TF ChIP-seq peak BED files (`best_site_sorted_unique.bed`, from `nucleosome_occupancy/`)
- MNase-seq dyad bedGraphs split by chromosome (from `nucleosome_occupancy/`)
- In vivo and in vitro `_best_score_nsm_profile.tab` files (from `nucleosome_occupancy/`)
- TF-to-profile mapping files (`avgs/tf_lists/tf_profile_mapping_<cell_line>.txt`)
- Reference genome FASTA (`hg38.fa`)

**Outputs**
- Part 1: 200 binned BED files per TF
- Part 2: `final_best_site_sorted_unique_{1..200}.bed` per TF
- Part 3: `best_site_sorted_unique_{i}_147bp_four_types.txt` per TF, per bin (the ΔNPS inputs)
- Averaged plots: `<cell_line>_hg38_<profile>_average_nsm_occ_deltaNPS_profile.png`
- Individual plots: `plots/individual/<cell_line>/<profile>/<tf>.png`
- Averaged ΔNPS data: `<cell_line>_<profile>_average_deltaNPS_profile_smoothed.txt`

---

## Dependencies

- Python ≥ 3.9 with numpy, matplotlib
- Perl ≥ 5.26 (for `separate_nucleosomes.pl`)
- bedtools ≥ 2.31
- bed2fasta 5.5.4
- SLURM (multi-node array jobs required for Part 2)
- Spack (optional, for environment management)

See the top-level `README.md` for tested versions and installation.
