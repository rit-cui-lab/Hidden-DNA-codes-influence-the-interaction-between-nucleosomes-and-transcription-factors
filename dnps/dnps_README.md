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

## Outputs

### `binning_bed.py` (Part 1)
200 BED files per TF, named `<tf>_<cell_line>_best_site_sorted_unique_<i>.bed` for *i* = 1..200. Each file contains the TF binding sites translated to one 10 bp window at the bin's offset.

```
chr22	41402983	41402992
chr22	42090033	42090042
```

### `callpipeline2.sh` (Part 2)
200 BED files per TF, named `final_best_site_sorted_unique_<i>.bed`, containing only the dyads that fall within each bin's windows. Real examples (bins 90–110) at [`demo/input/dnps_intermediates/`](../demo/input/dnps_intermediates/).

```
chr22	41403960	41403969	chr22	41403965	41403966	1.0
```

### `generate_147bp_bed.py` (Part 3)
A BED file (`<input>_147bp.bed`) where each row spans 147 bp around the dyad.

```
chr22	10510018	10510165
chr22	10510139	10510286
```

### `bed2fasta` (Part 3, external tool)
A FASTA file with one 147 bp sequence per nucleosome.

```
>chr22:10510018-10510165
ACGTACGTACGT...
```

### `separate_nucleosomes.pl` (Part 3)
A single-line file with four tab-separated integers — counts of nucleosomes classified as type 1, 2, 3, and 4 respectively. ΔNPS = (type1 − type4) / total × 100.

```
89	49	38	87
```

### `makeavgdnpsprofile.py`
A PNG plot (`<cell_line>_<profile>_average_nsm_occ_deltaNPS_profile.png`) with three y-axes (in vivo, ΔNPS, in vitro) plus a smoothed ΔNPS data file (`<cell_line>_<profile>_average_deltaNPS_profile_smoothed.txt`, 200 floats, one per line).

### `makeindividualplots.py`
A per-TF PNG plot. Also appends to a `tfs_without_invitro.txt` log when in vitro data is missing.

### dnps demo
A diffable TSV (`CTCF_dnps_demo.tsv`) plus a reference plot (`CTCF_dnps_demo.png`). Real examples at [`demo/expected/CTCF_dnps_demo.tsv`](../demo/expected/CTCF_dnps_demo.tsv) and [`demo/expected/CTCF_dnps_demo.png`](../demo/expected/CTCF_dnps_demo.png).

```
bin	distance_bp	type1	type2	type3	type4	total	delta_NPS_pct
90	-110	115	77	52	119	363	-1.101928
91	-100	86	63	45	93	287	-2.439024
```

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

Demo runs Part 3 only on 21 pre-computed bins (90–110) and produces a diffable TSV plus a reference plot.

---

## Dependencies

- Python ≥ 3.9 with numpy, matplotlib
- Perl ≥ 5.26
- bedtools ≥ 2.31
- bed2fasta 5.5.4
- SLURM (multi-node array jobs required for Part 2)
- Spack (optional, for environment management)

See the top-level `README.md` for tested versions and installation.
