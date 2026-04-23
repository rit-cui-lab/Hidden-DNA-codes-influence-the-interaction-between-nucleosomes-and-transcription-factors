# Nucleosome Occupancy Scripts Summary

## Overview
Computes genome-wide nucleosome occupancy from dyad position counts using Gaussian kernel smoothing, normalizes and converts to BigWig format, then generates occupancy profiles around transcription factor binding sites. This is the main pipeline that produces the in vivo and in vitro occupancy inputs consumed by `dnps/`.

---

## Script List

### Python Scripts

| File | Description |
|------|-------------|
| **calculatensm_multithread.py** | Computes nucleosome occupancy from dyad positions using a Gaussian kernel (σ = 20 bp, window ±73 bp). Uses `multiprocessing` to process chromosome chunks in parallel with overlap regions for accurate edges. Outputs bedGraph. |
| **calculate_occupancy_avg_and_normalize.py** | Normalizes occupancy by calculating the mean and dividing every score by it. Two-pass algorithm over the bedGraph. |
| **normalize_bedgraph.py** | Duplicate of `calculate_occupancy_avg_and_normalize.py` retained for compatibility with alternate bedGraph formats. Prefer `calculate_occupancy_avg_and_normalize.py` for new work. |
| **plotnucocc.py** | Creates nucleosome occupancy profile plots. Parses deepTools `plotProfile` output, applies smoothing, symmetrizes, and saves high-resolution PNG. Cell line selectable via CLI. |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **splitdyadfiles.sh** | Splits a full-genome dyad count file into per-chromosome files (chr1–22, chrX, chrY) for parallel downstream processing. |
| **sc_calculatensm_multitread.sh** | Batch SLURM submission for `calculatensm_multithread.py` across per-chromosome files. |
| **slurmcall_concatnsm.sh** | Concatenates per-chromosome occupancy bedGraphs into a single genome-wide file and compresses intermediates. |
| **slurmcallnormalizemnasedata_all.sh** | Batch SLURM submission for `calculate_occupancy_avg_and_normalize.py` across multiple cell lines. |
| **convert_bw.sh** | Generates per-cell-line SLURM jobs that sort, deduplicate, and convert the normalized bedGraph to BigWig via `bedGraphToBigWig`. |
| **convert_bwinvitro.sh** | Single-job variant of `convert_bw.sh` for the in vitro MNase-seq dataset. |
| **sc_nucocctfs.sh** | Generates and submits SLURM jobs for TF-centered occupancy analysis. For each MEME motif: BED→FASTA, FIMO motif scanning, bedtools merge, deepTools `computeMatrix` + `plotProfile`. |
| **sc_invitronucocctfs.sh** | In vitro variant of `sc_nucocctfs.sh`. Uses pre-computed in vitro BigWig and skips the FIMO step (reuses in vivo `best_site_sorted_unique.bed`). |
| **sc_plotnucocc.sh** | Generates per-TF SLURM jobs that call `plotnucocc.py` on the deepTools `.tab` profile output. |

---

## Pipeline / Data Flow

```
Dyad count file (from preprocessing/count_dyads.py)
    ↓
splitdyadfiles.sh                       (split by chromosome)
    ↓
sc_calculatensm_multitread.sh           → calculatensm_multithread.py
    ↓                                    (Gaussian-smoothed occupancy per chr)
slurmcall_concatnsm.sh                  (concatenate to genome-wide bedGraph)
    ↓
slurmcallnormalizemnasedata_all.sh      → calculate_occupancy_avg_and_normalize.py
    ↓                                    (divide by mean)
convert_bw.sh                            (bedGraph → BigWig)
    ↓
Normalized BigWig → consumed by dnps/ and TF-centered analysis below

TF-centered analysis (uses the normalized BigWig above plus MEME files):
    ↓
sc_nucocctfs.sh                          (per-TF FIMO + computeMatrix + plotProfile)
    ↓
Per-TF profile .tab + .png files
    ↓
sc_plotnucocc.sh → plotnucocc.py         (publication-quality plots)
```

---

## Usage

### Typical run order

```bash
# 1. Split the dyad file by chromosome
bash splitdyadfiles.sh <cell_line>_dyads_final.txt

# 2. Calculate occupancy per chromosome (edit CELL_LINE at the top first)
bash sc_calculatensm_multitread.sh <dir_of_per_chr_files>

# 3. Concatenate chromosome outputs
sbatch slurmcall_concatnsm.sh

# 4. Normalize
bash slurmcallnormalizemnasedata_all.sh

# 5. Convert to BigWig
bash convert_bw.sh

# 6. TF-centered analysis
bash sc_nucocctfs.sh <meme_dir> <bed_dir> <cell_line>
bash sc_invitronucocctfs.sh <meme_dir> <bed_dir> in_vitro

# 7. Final plots
bash sc_plotnucocc.sh
```

Edit the `CELL_LINE` variable, `#SBATCH` account/partition/email lines, and any hardcoded paths (e.g., `hg38.chrom.sizes`, `bedGraphToBigWig` location) in each script before submitting.

---

## Inputs and Outputs

**Inputs**
- Dyad count files (`<cell_line>_dyads_final.txt`) from `preprocessing/`
- Per-TF MEME files from `preprocessing/`
- ChIP-seq BED files (`.bed.gz`) from ENCODE
- Reference genome (`hg38.fa`) and chromosome sizes (`hg38.chrom.sizes`)
- In vitro MNase-seq BigWig (for in vitro analysis)

**Outputs**
- Per-chromosome smoothed occupancy bedGraphs
- Genome-wide concatenated, normalized bedGraph: `<cell_line>_nuc_occ_normalized.bg`
- BigWig version: `<cell_line>.bw`
- Per-TF directories under `chipseq/nucleosome_occupancy/<cell_line>/<tf>_<db>/` containing:
  - `best_site_sorted_unique.bed` — merged unique motif hit coordinates
  - `<tf>_<db>_best_site_score.gz` — deepTools matrix
  - `<tf>_<db>_best_score_nsm_profile.tab` — profile data (consumed by `dnps/` and `plotnucocc.py`)
  - `<cell_line>_hg38_<tf>_<db>_best_score_nsm_profile.png` — profile plot

---

## Dependencies

- Python ≥ 3.9 with numpy
- samtools ≥ 1.19
- bedtools ≥ 2.31
- bed2fasta 5.5.4
- deepTools (`computeMatrix`, `plotProfile`)
- MEME Suite (`fimo`)
- UCSC `bedGraphToBigWig` (vendored in this folder)
- SLURM (for batch submission)
- Spack (optional, for environment management)

See the top-level `README.md` for tested versions and installation.
