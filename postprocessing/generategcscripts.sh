#!/bin/bash
# generategcscripts.sh
# Author: Chris Carson
# Description: Generates batch SLURM job scripts for GC content analysis.
# Creates 200 individual job scripts per protein for BED→FASTA→GC pipeline.

# Configuration: Set cell line and base path
CELL_LINE=""
CELL_LINE_PATH="./chipseq/dNPS/${CELL_LINE}"

# Check if path exists
if [ ! -d "$CELL_LINE_PATH" ]; then
    echo "Error: Directory '$CELL_LINE_PATH' does not exist!"
    exit 1
fi

PROTEINS_TO_PROCESS=(
    ""
)

echo "Looking for proteins: ${PROTEINS_TO_PROCESS[@]}"
echo ""

# Reference to genome file (update path as needed)
GENOME_FILE="./genomes/hg38.fa"

for PROTEIN in "${PROTEINS_TO_PROCESS[@]}"; do
    
    # Find folder that starts with this protein name
    protein_dir=$(find "$CELL_LINE_PATH" -maxdepth 1 -type d -name "${PROTEIN}_*" | head -n 1)
    
    if [ -z "$protein_dir" ]; then
        echo "Warning: No folder found for protein $PROTEIN"
        continue
    fi
    
    # Add trailing slash
    protein_dir="${protein_dir}/"
    
    # Create job_scripts directory
    job_scripts_dir="${protein_dir}job_scripts/"
    mkdir -p "$job_scripts_dir"
    
    FOLDER_NAME=$(basename "$protein_dir")
    echo "Generating jobs for protein: $PROTEIN (from folder: $FOLDER_NAME)"
    
    # Generate 200 scripts for this protein
    for i in $(seq 1 200); do
        OUTPUT="${job_scripts_dir}${PROTEIN}_GC_${i}.sh"
        JOB_ID="${PROTEIN}_GC_${i}"
        
        cat > "$OUTPUT" << EOF
#!/bin/bash -l
#SBATCH --job-name=${JOB_ID}
#SBATCH -t 04:00:00
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15g
#SBATCH -o ${job_scripts_dir}/${JOB_ID}_output.txt
#SBATCH -e ${job_scripts_dir}/${JOB_ID}_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

spack env activate factor-x86_64-24062401

# Convert the dyad position to 10bp window
python ./scripts/postprocessing/generate_10bp_bed.py ${protein_dir}final_best_site_sorted_unique_${i}.bed

# Convert BED to FASTA
bed2fasta -s -both -o ${protein_dir}${PROTEIN}_best_site_sorted_unique_${i}_GC.bed.fa ${protein_dir}final_best_site_sorted_unique_${i}_10bp.bed ${GENOME_FILE}

# Calculate GC content of the sequences
python ./scripts/postprocessing/calculate_GC_content.py ${protein_dir}${PROTEIN}_best_site_sorted_unique_${i}_GC.bed.fa ${protein_dir}${PROTEIN}_best_site_sorted_unique_${i}_GC.txt
EOF

        chmod +x "$OUTPUT"
    done
    
    echo "  Generated 200 scripts in ${job_scripts_dir}"
    echo ""
done

echo "All done!"
