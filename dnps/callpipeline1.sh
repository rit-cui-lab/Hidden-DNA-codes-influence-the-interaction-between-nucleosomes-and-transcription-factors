#!/bin/bash
# =============================================================================
# callpipeline1.sh
# Author: Chris Carson
#
# Summary:
# This script generates and submits SLURM jobs for Part 1 of the dNPS pipeline.
# For each transcription factor MEME file, it creates 200 binned BED files that
# span -1000bp to +1000bp around TF binding site centers (10bp per bin).
# This binning is the foundation for calculating nucleosome positioning metrics
# at different distances from TF binding sites.
#
# Prerequisites:
#   - binning_bed.py script
#   - Input BED files with TF binding sites (best_site_sorted_unique.bed)
#   - Spack environment with required dependencies
#
# Output:
#   - 200 BED files per TF, each containing 10bp windows at a specific offset
# =============================================================================

# =============================================================================
# CONFIGURATION - Modify these paths for your environment
# =============================================================================
BASE_DIR="/path/to/your/project"
DNPS_DIR="${BASE_DIR}/chipseq/dNPS"
NUCOCC_DIR="${BASE_DIR}/chipseq/nucleosome_occupancy"
SCRIPTS_DIR="${BASE_DIR}/scripts/dnps"
JOB_SCRIPTS_DIR="${SCRIPTS_DIR}/job_scripts"
GENOME_FA="${BASE_DIR}/genomes/hg38.fa"
EMAIL="your_email@example.com"
ACCOUNT="your_account"
PARTITION="your_partition"
SPACK_ENV="your_spack_env"
# =============================================================================

# Create job scripts directory if it doesn't exist
mkdir -p "$JOB_SCRIPTS_DIR"

# Array of cell lines to process
CELL_LINES=("")

# Loop through each cell line
for cell_line in "${CELL_LINES[@]}"; do
    echo "Processing cell line: $cell_line"
    
    # Loop through each .meme directory in the cell line folder
    for meme_dir in "$NUCOCC_DIR/$cell_line"/*HOCO*.meme/; do
        [[ ! -d "$meme_dir" ]] && break
        
        meme_name=$(basename "$meme_dir")
        echo "  Processing meme directory: $meme_name"
        
        # Path to the bed file
        bed_file="${meme_dir}best_site_sorted_unique.bed"
        
        # Check if the bed file exists and create job script
        if [[ -f "$bed_file" ]]; then
            echo "    Creating job script for: $bed_file"
            
            # Create unique job ID and script name
            unique_id="${cell_line}_${meme_name}"
            job_script="${JOB_SCRIPTS_DIR}/${unique_id}.sh"
            output_dir="${DNPS_DIR}/${cell_line}/${meme_name}"
            mkdir -p "${output_dir}"

            # Generate the SLURM job script
            cat > "$job_script" << EOF
#!/bin/bash -l
#SBATCH -J ${unique_id}
#SBATCH -t 02:00:00
#SBATCH --mail-user=${EMAIL}
#SBATCH -A ${ACCOUNT} -p ${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15g
#SBATCH -o ${JOB_SCRIPTS_DIR}/${unique_id}_output.txt
#SBATCH -e ${JOB_SCRIPTS_DIR}/${unique_id}_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

spack env activate ${SPACK_ENV}

# Change to the meme directory
cd "$meme_dir"

# Create 200 bins for interval from -1000bp to +1000bp around TF motif, 10bp per bin
python ${SCRIPTS_DIR}/binning_bed.py best_site_sorted_unique.bed $meme_name $cell_line $output_dir
EOF

            # Make the job script executable
            chmod +x "$job_script"
            
            # Submit the job
            echo "    Submitting job: $unique_id"
            sbatch "$job_script"
            
            echo "    Job submitted for $bed_file"
        fi
    done
    
    echo "Completed processing cell line: $cell_line"
    echo "---"
done

echo "All jobs submitted!"
