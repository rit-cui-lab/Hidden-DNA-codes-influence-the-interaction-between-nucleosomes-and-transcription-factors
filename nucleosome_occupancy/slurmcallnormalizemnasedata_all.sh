#!/bin/bash
# slurmcallnormalizemnasedata_all.sh
# Author: Chris Carson
# Description: Batch SLURM job submission for normalizing MNase occupancy data.
# Normalizes occupancy across multiple cell lines by calculating and dividing by average.

# Configuration: Set cell lines and paths
CELL_LINES=("H1_ESC" "HepG2" "IMR90")
BASE_DIR="./mnase"
SCRIPT_PATH="./scripts/nucleosome_occupancy/calculate_occupancy_avg_and_normalize.py"

# Create directory for job scripts
JOB_SCRIPT_DIR="normalize_jobs"
mkdir -p "$JOB_SCRIPT_DIR"

for CELL_LINE in "${CELL_LINES[@]}"; do
    JOB_SCRIPT="${JOB_SCRIPT_DIR}/${CELL_LINE}_normalize_job.sh"
    
    cat > "$JOB_SCRIPT" << EOF
#!/bin/bash -l
#SBATCH -J ${CELL_LINE}-normalize
#SBATCH -t 2-7:0:0
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=25g
#SBATCH -o ${JOB_SCRIPT_DIR}/${CELL_LINE}_output.txt
#SBATCH -e ${JOB_SCRIPT_DIR}/${CELL_LINE}_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd ${BASE_DIR}/${CELL_LINE}

INPUT_FILE="${CELL_LINE}.bg"
OUTPUT_FILE="${CELL_LINE}_nuc_occ_normalized.bg"

echo "Normalizing occupancy for ${CELL_LINE}..."
python ${SCRIPT_PATH} \${INPUT_FILE} \${OUTPUT_FILE}

if [ \$? -eq 0 ]; then
    echo "Successfully normalized ${CELL_LINE}"
else
    echo "Error normalizing ${CELL_LINE}"
    exit 1
fi
EOF

    chmod +x "$JOB_SCRIPT"
    
    # Submit the job
    sbatch "$JOB_SCRIPT"
    echo "Submitted job for ${CELL_LINE}: $JOB_SCRIPT"
done

echo "All jobs submitted!"
echo "Job scripts are in: $JOB_SCRIPT_DIR"
