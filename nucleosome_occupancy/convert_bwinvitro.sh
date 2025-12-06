#!/bin/bash -l
# convert_bwinvitro.sh
# Author: Chris Carson
# Description: SLURM job for converting in vitro nucleosome occupancy to BigWig format.

#SBATCH -J convert-invitro
#SBATCH -t 07:00:00
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=25g
#SBATCH -o invitro_output.txt
#SBATCH -e invitro_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Configuration: Set paths relative to working directory
WORK_DIR="./mnase/in_vitro"
CHROM_SIZES="./scripts/nucleosome_occupancy/hg38.chrom.sizes"
BED_TO_BW_TOOL="./scripts/nucleosome_occupancy/bedGraphToBigWig"

cd "$WORK_DIR"

# Sort and deduplicate bedGraph
sort -k1,1 -k2,2n human_inVitro_sorted_dyad_normalized.bg | uniq > sorted_unique.bg

# Convert to BigWig
${BED_TO_BW_TOOL} sorted_unique.bg ${CHROM_SIZES} human_inVitro_sorted_dyad_normalized.bw

