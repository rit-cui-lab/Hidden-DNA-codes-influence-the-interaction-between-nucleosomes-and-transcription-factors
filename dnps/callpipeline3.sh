#!/bin/bash
# =============================================================================
# callpipeline3.sh
# Author: Chris Carson
#
# Summary:
# This script generates and submits SLURM jobs for Part 3 of the dNPS pipeline.
# For each of the 200 bins, it:
#   1. Converts nucleosome dyad positions to 147bp footprints
#   2. Extracts FASTA sequences for nucleosome regions
#   3. Classifies nucleosomes into four types based on dinucleotide patterns
#
# The four nucleosome types are used to calculate delta-NPS (Nucleosome 
# Positioning Score), which measures the directional preference of nucleosomes
# relative to TF binding sites.
#
# Prerequisites:
#   - Part 2 output (final_best_site_sorted_unique_{1-200}.bed files)
#   - Reference genome FASTA file
#   - bed2fasta tool
#   - separate_nucleosomes.pl script
#
# Output:
#   - 147bp BED files and FASTA sequences
#   - Four-types classification files for delta-NPS calculation
# =============================================================================

# =============================================================================
# CONFIGURATION - Modify these paths for your environment
# =============================================================================
BASE_DIR="/path/to/your/project"
DNPS_DIR="${BASE_DIR}/chipseq/dNPS"
SCRIPTS_DIR="${BASE_DIR}/scripts/dnps"
JOB_SCRIPTS_DIR="${SCRIPTS_DIR}/pt3_job_scripts"
GENOME_FA="${BASE_DIR}/genomes/hg38.fa"
EMAIL="your_email@example.com"
ACCOUNT="your_account"
PARTITION="your_partition"
SPACK_ENV="your_spack_env"
# =============================================================================

# Create job scripts directory if it doesn't exist
mkdir -p "$JOB_SCRIPTS_DIR"

# Array of cell lines
CELL_LINES=("HepG2")

# Loop through each cell line
for cell_line in "${CELL_LINES[@]}"; do
    echo "Processing cell line: $cell_line"

    # Loop through each .meme directory in the cell line folder
    for meme_dir in "$DNPS_DIR/$cell_line"/*.meme/; do
        [[ ! -d "$meme_dir" ]] && break

        meme_name=$(basename "$meme_dir")
        echo "  Processing meme directory: $meme_name"
        output_dir="${DNPS_DIR}/${cell_line}/${meme_name}"
        unique_id="${cell_line}_${meme_name}"
        job_script="${JOB_SCRIPTS_DIR}/part3_${unique_id}.sh"

        # Generate the SLURM job script
        cat > "$job_script" << EOF
#!/bin/bash -l
#SBATCH -J ${unique_id}
#SBATCH -t 03:00:00
#SBATCH --mail-user=${EMAIL}
#SBATCH -A ${ACCOUNT} -p ${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50g
#SBATCH -o ${JOB_SCRIPTS_DIR}/${unique_id}_output.txt
#SBATCH -e ${JOB_SCRIPTS_DIR}/${unique_id}_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

spack env activate ${SPACK_ENV}

# Change to the output directory
cd "${output_dir}"

for i in {1..200}; do
	# Convert dyad positions to 147bp nucleosome footprints
	python ${SCRIPTS_DIR}/generate_147bp_bed.py final_best_site_sorted_unique_\${i}.bed

	# Convert bed to fasta
	bed2fasta -s -both -o best_site_sorted_unique_\${i}_147bp.bed.fa final_best_site_sorted_unique_\${i}_147bp.bed ${GENOME_FA}

	# Calculate the number of nucleosomes with type1, type2, type3 and type4 patterns
	perl ${SCRIPTS_DIR}/separate_nucleosomes.pl best_site_sorted_unique_\${i}_147bp.bed.fa best_site_sorted_unique_\${i}_147bp_four_types.txt

done

EOF
        sbatch $job_script
    done
done
