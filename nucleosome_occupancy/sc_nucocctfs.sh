#!/bin/bash
# sc_nucocctfs.sh
# Author: Chris Carson
# Description: Generates SLURM jobs for nucleosome occupancy TF analysis.
# Matches MEME motifs with data files, runs FIMO, and generates profiles.

if [ $# -ne 3 ]; then
    echo "Usage: $0 <meme_directory> <data_directory> <cell_line>"
    exit 1
fi

MEME_DIR="$1"
DATA_DIR="$2"
CELL_LINE="$3"
BASE_DIR="./chipseq"
MNASE_DIR="./mnase"
GENOME_FILE="./genomes/hg38.fa"

# Create directories
JOB_SCRIPTS_DIR="${CELL_LINE}_nuc_occ_job_scripts"
mkdir -p "$JOB_SCRIPTS_DIR"

NUC_OCC_DIR="${BASE_DIR}/nucleosome_occupancy"
CELL_NUC_DIR="${NUC_OCC_DIR}/${CELL_LINE}"
mkdir -p "$CELL_NUC_DIR"

# Process MEME files
for meme_file in "$MEME_DIR"/*.meme; do
    echo "Processing MEME file: $meme_file"
    
    basename=$(basename "$meme_file")
    tf_name=${basename%%_*}
    db_name=$(echo "$basename" | cut -d'_' -f3)
    unique_id="${tf_name}_${db_name}"
    
    echo "Looking for bed files matching: $DATA_DIR/${tf_name}-*.bed.gz"
    
    # Match BED files
    for bed_file in "$DATA_DIR"/${tf_name}-*.bed.gz; do
        if [ -f "$bed_file" ]; then
            echo "Found matching bed file: $bed_file"
            
            tf_output_dir="${CELL_NUC_DIR}/${unique_id}"
            mkdir -p "$tf_output_dir"
            
            job_script="${JOB_SCRIPTS_DIR}/${unique_id}_job.sh"
            
            cat > "$job_script" << 'EOF'
#!/bin/bash -l
#SBATCH -J ${UNIQUE_ID}
#SBATCH -t 04:00:00
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
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
bw_file="${MNASE_DIR}/${CELL_LINE}/${CELL_LINE}.bw"
output_dir=${OUTPUT_DIR}
genome_file=${GENOME_FILE}

if [ ! -f "${bw_file}" ]; then
    echo "Error: Nucleosome occupancy file not found: ${bw_file}"
    exit 1
fi

cd ${output_dir}

# Create temp directory for BED file
temp_dir=$(mktemp -d)
gunzip -c ${bed_file} > ${temp_dir}/${unique_id}.bed

# Convert BED to FASTA
bed2fasta -s -both -o ${unique_id}.bed.fa ${temp_dir}/${unique_id}.bed ${genome_file}

# Create output directory for FIMO
mkdir -p fimo_output

echo "Running FIMO..."
fimo --oc fimo_output --verbosity 1 --bgfile --nrdb-- --text --thresh 1.0E-4 ${memefile} ${unique_id}.bed.fa > ${unique_id}_best_site.narrowPeak

# Rearrange and sort
awk 'BEGIN {OFS="\t"} {print $2, $3, $4, $1, $6, $5}' ${unique_id}_best_site.narrowPeak | tail -n +2 | bedtools sort -i > best_site_sorted.bed

# Merge overlapping regions
bedtools merge -i best_site_sorted.bed > best_site_sorted_unique.bed

echo "Computing matrix..."
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 -a 1000 \
    -R best_site_sorted_unique.bed \
    -S ${bw_file} \
    --skipZeros \
    -o ${unique_id}_best_site_score.gz \
    -p 6 \
    --outFileSortedRegions regions_${unique_id}_best_site_score.bed

echo "Plotting..."
plotProfile \
    -m ${unique_id}_best_site_score.gz \
    -out ${cell_line}_hg38_${unique_id}_best_score_nsm_profile.png \
    --outFileNameData ${unique_id}_best_score_nsm_profile.tab \
    --perGroup \
    --colors black \
    --plotTitle ${tf}_${db}_${cell_line} \
    --samplesLabel ${cell_line} \
    --refPointLabel "center" \
    -T ${tf}_${db} \
    -z ""

# Cleanup
rm -rf fimo_output ${temp_dir}
EOF

            # Substitute variables
            sed -i "s|\${UNIQUE_ID}|$unique_id|g" "$job_script"
            sed -i "s|\${TF_NAME}|$tf_name|g" "$job_script"
            sed -i "s|\${DB_NAME}|$db_name|g" "$job_script"
            sed -i "s|\${CELL_LINE}|$CELL_LINE|g" "$job_script"
            sed -i "s|\${MEME_FILE}|$meme_file|g" "$job_script"
            sed -i "s|\${BED_FILE}|$bed_file|g" "$job_script"
            sed -i "s|\${OUTPUT_DIR}|$tf_output_dir|g" "$job_script"
            sed -i "s|\${MNASE_DIR}|$MNASE_DIR|g" "$job_script"
            sed -i "s|\${GENOME_FILE}|$GENOME_FILE|g" "$job_script"
            sed -i "s|\${JOB_SCRIPTS_DIR}|$JOB_SCRIPTS_DIR|g" "$job_script"
            
            chmod +x "$job_script"
            sbatch "$job_script"
            echo "Created job script for ${unique_id}"
            break
        fi
    done
done

echo "Job script creation complete!"
