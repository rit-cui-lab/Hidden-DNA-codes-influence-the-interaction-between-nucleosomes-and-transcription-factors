#!/bin/bash
# convert_bw.sh
# Author: Chris Carson
# Description: Generates SLURM jobs for converting bedGraph to BigWig format.
# Creates one job per cell line with sorting and deduplication.

# Configuration: Set base directory and cell lines
BASE_DIR="./mnase"
CELL_LINES=(" ")
CHROM_SIZES="./scripts/nucleosome_occupancy/hg38.chrom.sizes"
BED_TO_BW_TOOL="./scripts/nucleosome_occupancy/bedGraphToBigWig"

# Create output directory for job scripts
JOB_SCRIPTS_DIR="${BASE_DIR}/convert_bw_scripts"
mkdir -p "$JOB_SCRIPTS_DIR"

for CELL_LINE in "${CELL_LINES[@]}"; do
    SCRIPT_NAME="${CELL_LINE}.sh"
    
    cat > "${SCRIPT_NAME}" << EOF
#!/bin/bash -l
#SBATCH -J ${CELL_LINE}
#SBATCH -t 07:00:00
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=25g
#SBATCH -o ${JOB_SCRIPTS_DIR}/${CELL_LINE}_output.txt
#SBATCH -e ${JOB_SCRIPTS_DIR}/${CELL_LINE}_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd ${BASE_DIR}/${CELL_LINE}

# Sort and deduplicate bedGraph file
sort -k1,1 -k2,2n ${CELL_LINE}_nuc_occ_normalized.bg | uniq > sorted_unique.bg

# Convert bedGraph to BigWig
${BED_TO_BW_TOOL} sorted_unique.bg ${CHROM_SIZES} ${CELL_LINE}.bw

EOF

    # Make executable and move to job scripts directory
    chmod +x "${SCRIPT_NAME}"
    mv "${SCRIPT_NAME}" "${JOB_SCRIPTS_DIR}/"
done

echo "Job script creation complete!"
