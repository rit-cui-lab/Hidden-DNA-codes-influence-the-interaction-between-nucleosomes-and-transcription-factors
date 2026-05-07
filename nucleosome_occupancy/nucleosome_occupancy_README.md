# Nucleosome Occupancy Scripts Summary

## Overview
Computes genome-wide nucleosome occupancy from dyad position counts using Gaussian kernel smoothing, normalizes and converts to BigWig format, then generates occupancy profiles around transcription factor binding sites. This is the main pipeline that produces the in vivo and in vitro occupancy inputs consumed by `dnps/`.

---

## Script List

### Python Scripts

| File | Description |
|------|-------------|
| **calculatensm_multithread.py** | Computes nucleosome occupancy from dyad positions using a Gaussian kernel (σ = 20 bp, window ±73 bp). Outputs bedGraph. |
| **calculate_occupancy_avg_and_normalize.py** | Normalizes occupancy by calculating the mean and dividing every score by it. |
| **normalize_bedgraph.py** | Duplicate of `calculate_occupancy_avg_and_normalize.py` retained for compatibility. Prefer `calculate_occupancy_avg_and_normalize.py` for new work. |
| **plotnucocc.py** | Creates nucleosome occupancy profile plots from deepTools `plotProfile` output. |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **splitdyadfiles.sh** | Splits a full-genome dyad count file into per-chromosome files. |
| **sc_calculatensm_multitread.sh** | Batch SLURM submission for `calculatensm_multithread.py`. |
| **slurmcall_concatnsm.sh** | Concatenates per-chromosome occupancy bedGraphs into one genome-wide file. |
| **slurmcallnormalizemnasedata_all.sh** | Batch SLURM submission for `calculate_occupancy_avg_and_normalize.py` across multiple cell lines. |
| **convert_bw.sh** | Per-cell-line SLURM jobs that sort, deduplicate, and convert the normalized bedGraph to BigWig. |
| **convert_bwinvitro.sh** | In vitro variant of `convert_bw.sh`. |
| **sc_nucocctfs.sh** | Per-MEME-motif SLURM jobs for TF-centered occupancy. |
| **sc_invitronucocctfs.sh** | In vitro variant of `sc_nucocctfs.sh`. |
| **sc_plotnucocc.sh** | Per-TF SLURM jobs that call `plotnucocc.py` on the deepTools `.tab` profile output. |

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

TF-centered analysis:
    ↓
sc_nucocctfs.sh                          (per-TF FIMO + computeMatrix + plotProfile)
    ↓
sc_plotnucocc.sh → plotnucocc.py         (publication-quality plots)
```

---

## Inputs

### `splitdyadfiles.sh` / `calculatensm_multithread.py`
Three-column TSV: chromosome, position, count (output of `preprocessing/count_dyads.py`). Real example: [`demo/input/HepG2_dyads_chr22.txt`](../demo/input/HepG2_dyads_chr22.txt).

```
chr22	10510092	2
chr22	10510213	1
```

### `calculate_occupancy_avg_and_normalize.py` / `normalize_bedgraph.py`
Four-column bedGraph: chromosome, start, end, score. Real example: [`demo/input/HepG2_chr22.bg`](../demo/input/HepG2_chr22.bg).

```
chr22	10510091	10510092	1.0
chr22	10510212	10510213	1.0026025852527254
```

### `convert_bw.sh`
The normalized bedGraph from the previous step plus a `hg38.chrom.sizes` file. Calls UCSC `bedGraphToBigWig`.

### `sc_nucocctfs.sh` / `sc_invitronucocctfs.sh`
A directory of per-TF MEME files, a directory of matching ENCODE ChIP-seq `.bed.gz` files, and a cell line name. Real examples: [`demo/input/CTCF_HUMAN_HOCOMOCOv11.meme`](../demo/input/CTCF_HUMAN_HOCOMOCOv11.meme) and [`demo/input/CTCF_HepG2_chr22.bed`](../demo/input/CTCF_HepG2_chr22.bed).

```
chr22	41403749	41404213	.	1000	.	681.83507	-1.00000	5.16371	225
```

### `plotnucocc.py` / `sc_plotnucocc.sh`
A `.tab` profile file from deepTools `plotProfile --outFileNameData`. First two lines are deepTools metadata; data lines are sample name followed by per-bin occupancy values.

---

## Usage

### Typical run order

```bash
bash nucleosome_occupancy/splitdyadfiles.sh HepG2_dyads_final.txt
bash nucleosome_occupancy/sc_calculatensm_multitread.sh per_chr_dir/
sbatch nucleosome_occupancy/slurmcall_concatnsm.sh
bash nucleosome_occupancy/slurmcallnormalizemnasedata_all.sh
bash nucleosome_occupancy/convert_bw.sh
bash nucleosome_occupancy/sc_nucocctfs.sh <meme_dir> <bed_dir> <cell_line>
bash nucleosome_occupancy/sc_plotnucocc.sh
```

Edit the `CELL_LINE` variable, `#SBATCH` account/partition/email lines, and any hardcoded paths in each script before submitting.

### Demo

```bash
bash nucleosome_occupancy/run_demo.sh demo/input demo/output/nucleosome_occupancy
```

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
