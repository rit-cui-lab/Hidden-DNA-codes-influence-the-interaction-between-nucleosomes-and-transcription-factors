# Preprocessing Scripts Summary

## Overview
Processes raw MNase-seq and ChIP-seq data into the dyad count files and MEME motif files used by the downstream pipelines (`nucleosome_occupancy`, `DAC`, `dnps`). Handles SAM/BAM filtering, dyad extraction, dyad counting, per-TF MEME file generation, and dataset organization.

---

## Script List

### Python Scripts

| File | Description |
|------|-------------|
| **count_dyads.py** | Counts and aggregates dyad positions from a tab-separated input file. Groups by chromosome and dyad position, outputting frequency counts sorted numerically by dyad position. |
| **generate_147bp_bed.py** | Converts dyad positions to BED format with 147 bp nucleosome windows (±73–74 bp around the dyad center). Produces standardized nucleosome coordinate ranges. |
| **processfactorbook.py** | Processes MEME motif files using ENCSR-to-protein mappings from a TSV file. Filters motifs by target protein list and optionally by biosample type. Generates one MEME file per protein. |
| **processmemefiles.py** | Extracts and organizes MEME motifs by protein name. Creates individual MEME files for each target protein with standardized naming. Tracks found and missing proteins. |

### Bash / SLURM Scripts

| File | Description |
|------|-------------|
| **filtersamfiles.sh** | Filters SAM files by insert size (120–180 bp), converts to BAM, sorts, and indexes. SLURM job script. |
| **generate_dyad_files_samfiles.sh** | Extracts dyad positions from paired-end SAM files by computing dyad coordinates from insert size and alignment start. SLURM job script. |
| **generate_dyad_files_samfiles_singleended.sh** | Extracts dyad positions from single-ended SAM files, handling both forward (flag=0) and reverse (flag=16) reads. SLURM job script. |
| **slurmcall_count_dyads.sh** | SLURM wrapper for `count_dyads.py`. |
| **renamebedfiles.sh** | Renames ChIP-seq BED files based on experiment target metadata from a TSV file. |
| **removesharedmemes.sh** | Compares two folders of MEME files and removes from folder 2 any proteins already present in folder 1. Prompts for confirmation before deletion. |
| **sc_processmemefiles.sh** | Batch SLURM submission for `processmemefiles.py` across a directory of MEME files. |
| **compareproteins.sh** *(verify file exists in repo)* | Compares protein datasets between two directories of MEME files. Reports proteins shared, unique to folder 1, and unique to folder 2. |
| **slurmcall_processmemefiles_all.sh** *(verify file exists in repo)* | Batch SLURM submission that generates per-file job scripts for processing multiple MEME files. |

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

ChIP-seq BED files (from ENCODE)
    ↓
renamebedfiles.sh            (rename using metadata.tsv)
    ↓
Renamed BED files → consumed by nucleosome_occupancy/
```

---

## Usage

### Typical run order

```bash
# 1. Filter SAM and extract dyads (edit CELL_LINE inside each script first)
sbatch filtersamfiles.sh
sbatch generate_dyad_files_samfiles.sh        # or the single-ended variant

# 2. Count dyads
sbatch slurmcall_count_dyads.sh

# 3. Prepare MEME files for each TF
bash sc_processmemefiles.sh <meme_dir> <protein_list.txt> <cell_line>

# 4. Rename downloaded ChIP-seq BED files
bash renamebedfiles.sh
```

Edit the `CELL_LINE` variable and the `#SBATCH` account/partition/email lines at the top of each script before submitting.

---

## Inputs and Outputs

**Inputs**
- Aligned SAM files (paired-end or single-ended) — one per cell line
- MEME motif files from HOCOMOCO v11 and FactorBook
- ChIP-seq BED files with accompanying `metadata.tsv` from ENCODE
- Protein list text file (one gene symbol per line)

**Outputs**
- `filtered_<cell_line>.sam` — insert-size-filtered reads
- `<cell_line>_sorted.bam` + `.bai` — sorted indexed BAM
- `<cell_line>_dyads.txt` — raw dyad positions
- `<cell_line>_dyads_final.txt` — dyad positions with counts (chromosome, position, count)
- `<protein>_HUMAN_<source>.meme` — one MEME file per target protein
- Renamed `<TF>-<file_id>.bed.gz` ChIP-seq BED files

---

## Dependencies

- Python ≥ 3.9 (standard library only — no pip packages)
- Perl ≥ 5.26
- samtools ≥ 1.19
- SLURM (for batch submission)
- Spack (optional, for environment management on HPC)

See the top-level `README.md` for tested versions and installation.
