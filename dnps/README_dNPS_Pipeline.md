# dNPS Pipeline Scripts
## Author: Chris Carson

This collection of scripts implements the delta-NPS (Nucleosome Positioning Score) pipeline for analyzing nucleosome positioning around transcription factor binding sites.

---

## File Summary

### Python Scripts

| File | Description |
|------|-------------|
| **binning_bed.py** | Creates 200 binned BED files from TF binding sites. Each bin represents a 10bp window spanning -1000bp to +1000bp around the binding site center. This is the foundation for spatial analysis in Part 1 of the pipeline. |
| **generate_147bp_bed.py** | Converts nucleosome dyad positions to 147bp nucleosome footprints. Takes single-coordinate dyad positions and expands them to full nucleosome coverage (-74bp to +73bp). Used in Part 3 of the pipeline. |
| **makeavgdnpsprofile.py** | Generates averaged nucleosome positioning profile plots for groups of TFs. Combines in vivo occupancy (red), delta-NPS (blue), and in vitro occupancy (green) data across multiple TFs of the same profile category. |
| **makeindividualplots.py** | Generates individual nucleosome positioning profile plots for single TFs. Similar to the averaging script but produces per-TF visualizations. Handles missing in vitro data gracefully. |

### Bash/SLURM Scripts

| File | Description |
|------|-------------|
| **callpipeline1.sh** | Part 1 of dNPS pipeline. Submits SLURM jobs to run binning_bed.py for each TF, creating the 200 binned BED files needed for subsequent analysis. |
| **callpipeline2.sh** | Part 2 of dNPS pipeline. Submits multi-node SLURM array jobs to intersect binned TF regions with MNase-seq nucleosome dyad positions using bedtools. |
| **callpipeline3.sh** | Part 3 of dNPS pipeline. Submits SLURM jobs to convert dyad positions to 147bp footprints, extract FASTA sequences, and classify nucleosomes into four types for delta-NPS calculation. |
| **makeavgplotscripts.sh** | Submits SLURM jobs to create averaged profile plots. Organizes TFs by profile category (DIP/DIP, DIP/PEAK, PEAK/DIP, PEAK, AMBIGUOUS) and generates group-averaged visualizations. |
| **callindividualplots.sh** | Submits SLURM jobs to create individual profile plots for each TF. Plots are organized by cell line and profile category. |
| **call3lineplots.sh** | Alternative plotting script submission. Similar to makeavgplotscripts.sh but calls make3lineplots.py for potentially different visualization options. |

---

## Pipeline Overview

```
Part 1 (callpipeline1.sh)
    │
    ├── binning_bed.py
    │   └── Creates 200 binned BED files per TF
    │
    ▼
Part 2 (callpipeline2.sh)
    │
    ├── bedtools intersect
    │   └── Intersects bins with MNase-seq dyad positions
    │
    ▼
Part 3 (callpipeline3.sh)
    │
    ├── generate_147bp_bed.py
    │   └── Converts dyads to 147bp footprints
    │
    ├── bed2fasta
    │   └── Extracts nucleosome sequences
    │
    ├── separate_nucleosomes.pl
    │   └── Classifies nucleosomes into 4 types
    │
    ▼
Visualization
    │
    ├── makeavgplotscripts.sh → makeavgdnpsprofile.py
    │   └── Averaged plots by profile category
    │
    ├── callindividualplots.sh → makeindividualplots.py
    │   └── Individual TF plots
    │
    └── call3lineplots.sh → make3lineplots.py
        └── Alternative averaged plots
```

---

## Configuration

All scripts contain a **CONFIGURATION** section at the top where you should modify paths for your environment:

```bash
# =============================================================================
# CONFIGURATION - Modify these paths for your environment
# =============================================================================
BASE_DIR="/path/to/your/project"
EMAIL="your_email@example.com"
ACCOUNT="your_account"
PARTITION="your_partition"
SPACK_ENV="your_spack_env"
# =============================================================================
```

---

## Profile Categories

The pipeline categorizes TFs into five profile types based on their nucleosome positioning patterns:

- **DIP_DIP** - Dip in both in vivo and in vitro occupancy
- **DIP_PEAK** - Dip in vivo, peak in vitro
- **PEAK_DIP** - Peak in vivo, dip in vitro
- **PEAK** - Peak in both datasets
- **AMBIGUOUS** - Unclear or mixed patterns

---

## Dependencies

- Python 3.x with numpy, matplotlib
- bedtools
- bed2fasta
- Perl (for separate_nucleosomes.pl)
- SLURM workload manager
- Spack package manager (optional, for environment management)

---

## Output

The pipeline produces:
1. Binned BED files for each distance from TF center
2. Nucleosome classification files (four types)
3. PNG plot files showing nucleosome positioning profiles
