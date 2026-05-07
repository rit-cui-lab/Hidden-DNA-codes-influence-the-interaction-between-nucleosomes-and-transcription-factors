# Preprocessing Scripts Summary

## Overview
Processes raw MNase-seq and ChIP-seq data into the dyad count files and MEME motif files used by the downstream pipelines (`nucleosome_occupancy`, `DAC`, `dnps`). Handles SAM/BAM filtering, dyad extraction, dyad counting, per-TF MEME file generation, and dataset organization.

---

## Script List

### Python Scripts

| File | Description |
|------|-------------|
| **count_dyads.py** | Counts and aggregates dyad positions from a tab-separated input file. |
| **generate_147bp_bed.py** | Converts dyad positions to BED format with 147 bp nucleosome windows (±73–74 bp around the dyad center). |
| **processfactorbook.py** | Processes MEME motif files using ENCSR-to-protein mappings from a TSV file. Filters motifs by target protein list and optionally by biosample type. |
| **processmemefiles.py** | Extracts and organizes MEME motifs by protein name. Creates individual MEME files for each target protein. |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **filtersamfiles.sh** | Filters SAM files by insert size (120–180 bp), converts to BAM, sorts, and indexes. SLURM job script. |
| **generate_dyad_files_samfiles.sh** | Extracts dyad positions from paired-end SAM files. SLURM job script. |
| **generate_dyad_files_samfiles_singleended.sh** | Extracts dyad positions from single-ended SAM files. SLURM job script. |
| **slurmcall_count_dyads.sh** | SLURM wrapper for `count_dyads.py`. |
| **renamebedfiles.sh** | Renames ChIP-seq BED files based on experiment target metadata from a TSV file. |
| **removesharedmemes.sh** | Compares two folders of MEME files and removes from folder 2 any proteins already present in folder 1. |
| **sc_processmemefiles.sh** | Batch SLURM submission for `processmemefiles.py` across a directory of MEME files. |

---

## Pipeline / Data Flow

```
SAM files (from upstream bowtie2 alignment)
    ↓
filtersamfiles.sh           (insert-size filter, SAM→BAM, sort, index)
    ↓
generate_dyad_files_samfiles[_singleended].sh
    ↓                        (extract dyad positions per read)
count_dyads.py → slurmcall_count_dyads.sh
    ↓                        (aggregate duplicate dyads into counts)
Dyad count files → consumed by DAC/, nucleosome_occupancy/, dnps/

MEME motif files (from HOCOMOCO / FactorBook / ENCODE)
    ↓
processmemefiles.py / processfactorbook.py
    ↓                        (split into one MEME file per TF)
Per-TF MEME files → consumed by nucleosome_occupancy/, dnps/
```

---

## Inputs

### `filtersamfiles.sh` / `generate_dyad_files_samfiles.sh`
Standard paired-end SAM from bowtie2. Pipeline reads columns 2 (flag), 3 (chromosome), 4 (position), and 9 (insert size).

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr22	LN:50818468
read_001	99	chr22	10510055	60	147M	=	10510155	147	ACGT...	IIII...
read_001	147	chr22	10510155	60	147M	=	10510055	-147	ACGT...	IIII...
```

### `count_dyads.py`
Two-column TSV: chromosome and dyad position. One line per read. No header.

```
chr22	10510092
chr22	10510092
chr22	10510213
```

### `generate_147bp_bed.py`
Three-column TSV: chromosome, dyad position, count (output of `count_dyads.py`). Real example: [`demo/input/HepG2_dyads_chr22.txt`](../demo/input/HepG2_dyads_chr22.txt).

```
chr22	10510092	2
chr22	10510213	1
```

### `processmemefiles.py` / `processfactorbook.py`
Standard MEME motif file (concatenated motifs) plus a plain-text protein list (one gene symbol per line). `processfactorbook.py` additionally requires an ENCSR-to-protein TSV mapping. See [`demo/input/CTCF_HUMAN_HOCOMOCOv11.meme`](../demo/input/CTCF_HUMAN_HOCOMOCOv11.meme) for a single-TF MEME example.

### `renamebedfiles.sh`
A folder of `.bed.gz` files plus a `metadata.tsv` from ENCODE in the same folder. Script reads column 23 of the TSV (experiment target).

---

## Usage

### Typical run order

```bash
sbatch preprocessing/filtersamfiles.sh
sbatch preprocessing/generate_dyad_files_samfiles.sh
sbatch preprocessing/slurmcall_count_dyads.sh
bash preprocessing/sc_processmemefiles.sh <meme_dir> <protein_list.txt> <cell_line>
bash preprocessing/renamebedfiles.sh
```

Edit the `CELL_LINE` variable and the `#SBATCH` account/partition/email lines at the top of each script before submitting.

### Demo

```bash
bash preprocessing/run_demo.sh demo/input demo/output/preprocessing
```

---

## Dependencies

- Python ≥ 3.9 (standard library only — no pip packages)
- Perl ≥ 5.16
- samtools ≥ 1.19
- SLURM (for batch submission)
- Spack (optional, for environment management on HPC)

See the top-level `README.md` for tested versions and installation.
