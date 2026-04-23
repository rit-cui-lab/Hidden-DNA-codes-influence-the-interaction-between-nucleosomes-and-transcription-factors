# Hidden DNA codes influence the interaction between nucleosomes and transcription factors

Code and pipelines accompanying the manuscript *"Intrinsic DNA codes govern
distinct modes of nucleosome–transcription factor interactions"* (Carson et al.).

Preprint: <https://www.biorxiv.org/content/10.64898/2025.12.30.697010v1>

Run order: `preprocessing` → `nucleosome_occupancy` + `DAC` + `dnps` + `domain_analysis` → `postprocessing`.

## 1. System requirements

**Operating systems tested:**
- Red Hat Enterprise Linux 9 (RIT SPORC HPC cluster)

**Languages and interpreters:**
- Python 3.9 (tested with 3.9.25)
- R 4.5 (tested with 4.5.2)
- Perl ≥ 5.16 (for `DAC/auto_correlation_multiple_new4_corrected_new.pl` and `separate_nucleosomes.pl`)
- GNU bash ≥ 4.4

**Python packages** (see `requirements.txt`): numpy, matplotlib, pandas, requests, biopython.

**R packages** (see `install.R`): dplyr, ggplot2, tidyr, stringr, tibble.

**External command-line tools:**
- samtools 1.19.2 (htslib 1.19.1) — used in `preprocessing/filtersamfiles.sh` for SAM/BAM conversion, sorting, and indexing
- bedtools 2.31.0 — used in `dnps/callpipeline2.sh` for interval intersection
- bed2fasta 5.5.4 — used in `dnps/callpipeline3.sh` and `postprocessing/generategcscripts.sh` for sequence extraction
- deepTools (`computeMatrix`, `plotProfile`) — used in `nucleosome_occupancy/sc_nucocctfs.sh` and `sc_invitronucocctfs.sh` for occupancy profile generation
- MEME Suite (`fimo`) — used in `nucleosome_occupancy/sc_nucocctfs.sh` for motif scanning
- UCSC `bedGraphToBigWig` — used in `nucleosome_occupancy/convert_bw.sh` for bedGraph→BigWig conversion (vendored in `nucleosome_occupancy/`; source: <http://hgdownload.soe.ucsc.edu/admin/exe/>)

**Pipeline inputs:** This pipeline begins from aligned SAM files. Read alignment (MNase-seq → SAM) is performed upstream using bowtie2 and is not included in this repository; alignment parameters are described in the Methods section of the manuscript.

The original environment was built with Spack on an HPC cluster and jobs were submitted via Slurm. Reviewers do not need Spack or Slurm to run the demo — any installation of the tools above (pip, conda, system packages) is sufficient.

**Hardware:** ≥ 16 GB RAM recommended. No GPU required. A Slurm cluster is recommended for full-genome analysis but not for the demo.

## 2. Installation guide

```bash
git clone https://github.com/rit-cui-lab/Hidden-DNA-codes-influence-the-interaction-between-nucleosomes-and-transcription-factors.git
cd Hidden-DNA-codes-influence-the-interaction-between-nucleosomes-and-transcription-factors
pip install -r requirements.txt
Rscript install.R
```

External tools (`samtools`, `bedtools`, `bed2fasta`) must be installed separately and available on `PATH`. Note that `bed2fasta` is not in most standard bioinformatics distributions and may require separate sourcing; the `dnps` and `postprocessing` demos will skip gracefully if it is not present.

**Typical install time on a normal desktop computer:** ~2 minutes for the Python and R packages. External command-line tools install time varies by system package manager.

## 3. Demo

A small demo dataset is provided in `demo/input/` (HepG2 MNase-seq on chromosome 22, a subset of CTCF ChIP-seq peaks on chromosome 22, and a chromosome 22 reference FASTA). Each pipeline stage has its own `run_demo.sh` that runs without Slurm.

**Prerequisites for the dnps and postprocessing demos:** these two demos require `bed2fasta` (tested with 5.5.4) on `PATH`. `bed2fasta` is not a standard bioinformatics distribution tool and may not be available on a typical workstation. On the RIT SPORC cluster it is provided by `spack env activate factor-x86_64-24062401`. On other systems, install it separately and ensure it is on `PATH` before running those demos. The preprocessing, nucleosome_occupancy, and DAC demos have no such requirement and will run on any system with Python, Perl, numpy, and matplotlib installed.

After running, reviewers can verify reproducibility by diffing the outputs against the reference files in `demo/expected/`.

### preprocessing demo (~5 seconds)

```bash
bash preprocessing/run_demo.sh demo/input demo/output/preprocessing
```

Runs `count_dyads.py` on a chr22 slice of HepG2 dyad positions (expanded back to one line per occurrence). Output:
- `HepG2_dyads_chr22_counts.txt` — dyad positions with counts (compare against `demo/expected/HepG2_dyads_chr22_counts.head.txt`).

### nucleosome_occupancy demo (~1 minute)

```bash
bash nucleosome_occupancy/run_demo.sh demo/input demo/output/nucleosome_occupancy
```

Runs `calculatensm_multithread.py` on the first 50,000 chr22 dyads, then normalizes with `calculate_occupancy_avg_and_normalize.py`. Output:
- `HepG2_dyads_chr22_demo.bg` — raw Gaussian-smoothed occupancy
- `HepG2_dyads_chr22_demo_normalized.bg` — normalized by average (compare against `demo/expected/HepG2_chr22_smoothed.head.bg`, first 100 lines)

### DAC demo (~1 minute)

```bash
bash DAC/run_demo.sh demo/input demo/output/DAC
```

Runs the Perl distance auto-correlation on the first 20,000 chr22 dyads. Output:
- `HepG2_chr22_auto_correlation.txt` — distance vs auto-correlation count histogram (compare against `demo/expected/CTCF_HepG2_chr22_DAC.txt`).

### dnps demo (~30 seconds, requires bed2fasta)

```bash
bash dnps/run_demo.sh demo/input demo/output/dnps
```

Runs dnps Part 3 (147 bp nucleosome classification) on 21 pre-computed bins (bins 90–110, ±100 bp around CTCF binding site centers), then aggregates into a diffable table. Output:
- `CTCF_dnps_demo.tsv` — tab-separated table of type1–4 counts and ΔNPS per bin; this is the primary deliverable for reproducibility checks (`diff` against `demo/expected/CTCF_dnps_demo.tsv`)
- `CTCF_dnps_demo.png` — reference plot (unsmoothed; noisy at this sample size, for visual inspection only)

### postprocessing demo (~20 seconds, requires bed2fasta)

```bash
bash postprocessing/run_demo.sh demo/input demo/output/postprocessing
```

Computes GC content for bins 90–110 around CTCF centers and generates a plot. Output:
- `CTCF_chr22_GC.png` — smoothed GC content profile (compare against `demo/expected/CTCF_chr22_GC.png`)

### Expected total run time

~3 minutes on a normal desktop computer (no Slurm, no HPC resources required).

## 4. Instructions for use

To run on your own data:

1. **Inputs:** aligned BAM files (MNase-seq, ChIP-seq), a BED file of TF binding sites, and a reference genome FASTA.
2. Edit the `CONFIGURATION` blocks at the top of each `.sh` script and the `BASE_DIR` constants in `dnps/makeavgdnpsprofile.py` and `dnps/makeindividualplots.py`.
3. Run the pipeline directories in the order listed at the top of this README. On an HPC cluster with Slurm, the `call*.sh` scripts in each directory auto-generate and submit jobs — adjust the `#SBATCH` account/partition/email values for your cluster first.

## Citation

If you use this code, please cite:

> Carson, C. W., Umesh Nagalakshmi, S., Adhikari, I., Freewoman, J. M.,
> Pizzi, J. R., Prakash, P., & Cui, F. (2025). Intrinsic DNA codes govern
> distinct modes of nucleosome–transcription factor interactions. *bioRxiv*.
> <https://doi.org/10.64898/2025.12.30.697010>

## License

Released under the [MIT License](LICENSE).

## Contact

**Christopher Carson** (lead developer) — <christophercarsonscience@gmail.com> — University at Buffalo, Jacobs School of Medicine and Biomedical Sciences.

**Feng Cui** (PI) — <fxcsbi@rit.edu> — Thomas H. Gosnell School of Life Sciences, Rochester Institute of Technology.
