#!/bin/bash
# sc_invitronucocctfs.sh
# Author: Chris Carson
# Description: Generates SLURM jobs for in vitro nucleosome occupancy analysis.
# Matches MEME motif files with data files and generates occupancy profiles.

if [ $# -ne 3 ]; then
    echo "Usage: $0 <meme_directory> <data_directory> <cell_line>"
    echo "Example: $0 ./chipseq/memes ./chipseq/bed_files in_vitro"
    exit 1
fi

MEME_DIR="$1"
DATA_DIR="$2"
CELL_LINE="$3"

# Configuration: Set base directories
BASE_DIR="./chipseq"
MNASE_DIR="./mnase/in_vitro"
GENOME_FILE="./genomes/hg38.fa"

# Create directories for job scripts and outputs
JOB_SCRIPTS_DIR="${CELL_LINE}_nuc_occ_job_scripts"
mkdir -p "$JOB_SCRIPTS_DIR"

NUC_OCC_DIR="${BASE_DIR}/nucleosome_occupancy"
CELL_NUC_DIR="${NUC_OCC_DIR}/${CELL_LINE}"
mkdir -p "$CELL_NUC_DIR"

echo "Processing MEME files from: $MEME_DIR"
echo "Data files from: $DATA_DIR"
echo "Cell line: $CELL_LINE"
echo "Job scripts directory: $JOB_SCRIPTS_DIR"
echo "Output directory: $CELL_NUC_DIR"
echo ""

JOB_COUNT=0

# Process each MEME file
for meme_file in "$MEME_DIR"/*.meme; do
    # Skip if no files found
    [ -e "$meme_file" ] || continue
    
    echo "Processing MEME file: $(basename "$meme_file")"
    
    # Extract TF and database names from filename
    basename=$(basename "$meme_file")
    tf_name=${basename%%_*}
    db_name=$(echo "$basename" | cut -d'_' -f3)
    unique_id="${tf_name}_${db_name}"
    
    # Look for matching BED files
    for bed_file in "$DATA_DIR"/${tf_name}-*.bed.gz; do
        # Verify BED file exists
        if [ -f "$bed_file" ]; then
            echo "  Found matching BED file: $(basename "$bed_file")"
            
            # Create output directory for this TF
            tf_output_dir="${CELL_NUC_DIR}/${unique_id}"
            mkdir -p "$tf_output_dir"
            
            # Create SLURM job script
            job_script="${JOB_SCRIPTS_DIR}/${unique_id}_job.sh"
            
            cat > "$job_script" << 'EOF'
#!/bin/bash -l
#SBATCH -J ${UNIQUE_ID}
#SBATCH -t 04:00:00
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH -o ${JOB_SCRIPTS_DIR}/${UNIQUE_ID}_output.txt
#SBATCH -e ${JOB_SCRIPTS_DIR}/${UNIQUE_ID}_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

spack env activate factor-x86_64-24062401

tf=${TF_NAME}
db=${DB_NAME}
unique_id=${UNIQUE_ID}
cell_line=${CELL_LINE}
memefile=${MEME_FILE}
bed_file=${BED_FILE}
bw_file="${MNASE_DIR}/human_inVitro_sorted_dyad_normalized.bw"
output_dir=${OUTPUT_DIR}

# Verify nucleosome occupancy file exists
if [ ! -f "${bw_file}" ]; then
    echo "Error: Nucleosome occupancy file not found: ${bw_file}"
    exit 1
fi

cd ${output_dir}

echo "Computing matrix for ${tf}_${db}..."
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 -a 1000 \
    -R best_site_sorted_unique.bed \
    -S ${bw_file} \
    --skipZeros \
    -o invitro_${unique_id}_best_site_score.gz \
    -p 6 \
    --outFileSortedRegions invitro_regions_${unique_id}_best_site_score.bed

echo "Generating plot..."
plotProfile \
    -m invitro_${unique_id}_best_site_score.gz \
    -out invitro_${cell_line}_hg38_${unique_id}_best_score_nsm_profile.png \
    --outFileNameData invitro_${unique_id}_best_score_nsm_profile.tab \
    --perGroup \
    --colors black \
    --plotTitle ${tf}_${db}_${cell_line} \
    --samplesLabel ${cell_line} \
    --refPointLabel "center" \
    -T ${tf}_${db} \
    -z ""

echo "Complete!"
EOF

            # Substitute variables in the script
            sed -i "s|\${UNIQUE_ID}|$unique_id|g" "$job_script"
            sed -i "s|\${TF_NAME}|$tf_name|g" "$job_script"
            sed -i "s|\${DB_NAME}|$db_name|g" "$job_script"
            sed -i "s|\${CELL_LINE}|$CELL_LINE|g" "$job_script"
            sed -i "s|\${MEME_FILE}|$meme_file|g" "$job_script"
            sed -i "s|\${BED_FILE}|$bed_file|g" "$job_script"
            sed -i "s|\${OUTPUT_DIR}|$tf_output_dir|g" "$job_script"
            sed -i "s|\${JOB_SCRIPTS_DIR}|$JOB_SCRIPTS_DIR|g" "$job_script"
            
            chmod +x "$job_script"
            
            # Submit the job
            sbatch "$job_script"
            echo "  Created and submitted job: $job_script"
            ((JOB_COUNT++))
            
            # Only process first matching BED file
            break
        fi
    done
done

echo ""
echo "========================================="
echo "Job script creation complete!"
echo "Total jobs created: $JOB_COUNT"
echo "Job scripts directory: $JOB_SCRIPTS_DIR"
echo "Output directory: $CELL_NUC_DIR"
echo "========================================="
