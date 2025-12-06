#!/bin/bash
# =============================================================================
# makeavgplotscripts.sh
# Author: Chris Carson
#
# Summary:
# This script generates and submits SLURM jobs to create averaged nucleosome
# positioning profile plots. It organizes transcription factors by their
# profile categories (DIP/DIP, DIP/PEAK, PEAK/DIP, PEAK, AMBIGUOUS) and creates
# averaged plots for each category showing:
#   - In vivo nucleosome occupancy (red)
#   - Delta-NPS percentage (blue)
#   - In vitro nucleosome occupancy (green)
#
# The script reads TF-to-profile mappings from cell-line-specific files and
# groups TFs accordingly before generating the visualization jobs.
#
# Prerequisites:
#   - Completed dNPS pipeline (parts 1-3)
#   - Nucleosome occupancy profile files
#   - TF-to-profile mapping files
#   - makeavgdnpsprofile.py script
#
# Output:
#   - Averaged profile plots per cell line and profile category
# =============================================================================

# =============================================================================
# CONFIGURATION - Modify these paths for your environment
# =============================================================================
BASE_DIR="/path/to/your/project"
DNPS_DIR="${BASE_DIR}/chipseq/dNPS"
NUCOCC_DIR="${BASE_DIR}/chipseq/nucleosome_occupancy"
SCRIPTS_DIR="${BASE_DIR}/scripts/dnps"
JOB_SCRIPTS_DIR="${SCRIPTS_DIR}/plot_job_scripts/avgs"
TF_LISTS_DIR="${JOB_SCRIPTS_DIR}/tf_lists"
TF_PROFILE_MAPPING_DIR="${SCRIPTS_DIR}/avgs/tf_lists"
PLOTS_DIR="${DNPS_DIR}/plots"
EMAIL="your_email@example.com"
ACCOUNT="your_account"
PARTITION="your_partition"
# =============================================================================

# Create directories if they don't exist
mkdir -p "$JOB_SCRIPTS_DIR"
mkdir -p "$TF_LISTS_DIR"

# Array of cell lines
CELL_LINES=("MCF-7")

# Array of profiles
PROFILES=("DIP_DIP" "DIP_PEAK" "PEAK_DIP" "PEAK" "AMBIGUOUS")

# Loop through each cell line
for cell_line in "${CELL_LINES[@]}"; do
    echo "Processing cell line: $cell_line"
    
    # Use cell-line-specific mapping file
    tf_profile_mapping="${TF_PROFILE_MAPPING_DIR}/tf_profile_mapping_${cell_line}.txt"
    
    # Check if mapping file exists
    if [[ ! -f "$tf_profile_mapping" ]]; then
        echo "ERROR: Mapping file not found: $tf_profile_mapping"
        continue
    fi
    
    echo "Using mapping file: $tf_profile_mapping"
    
    # Initialize TF list files for each profile
    declare -A profile_files
    for profile in "${PROFILES[@]}"; do
        tf_list_file="${TF_LISTS_DIR}/${cell_line}_${profile}_tf_list.txt"
        profile_files[$profile]=$tf_list_file
        > "$tf_list_file"
    done
    
    # Collect all TFs for this cell line and categorize by profile
    for meme_dir in "$DNPS_DIR/$cell_line"/*.meme/; do
        [[ ! -d "$meme_dir" ]] && break
        
        meme_name=$(basename "$meme_dir" .meme)
        
        # Extract just the TF name (remove _HOCOMOCOv11 or _complete-factorbook-catalog suffix)
        tf_name=$(echo "$meme_name" | sed 's/_HOCOMOCOv11$//' | sed 's/_complete-factorbook-catalog$//')
        
        nucocc_file="${NUCOCC_DIR}/${cell_line}/${meme_name}.meme/${meme_name}.meme_best_score_nsm_profile.tab"
        
        # Check if nucleosome occupancy file exists
        if [[ -f "$nucocc_file" ]]; then
            # Determine which profile this TF belongs to using the extracted TF name
            profile=$(awk -v tf="${tf_name}" '$1 == tf {print $2}' "$tf_profile_mapping")
            
            # Strip any trailing whitespace/newlines
            profile=$(echo "$profile" | tr -d '\n\r')
            
            if [[ -n "$profile" ]] && [[ " ${PROFILES[@]} " =~ " ${profile} " ]]; then
                echo "DEBUG: About to write meme_name='$meme_name'" 
                echo "$nucocc_file $meme_name" >> "${profile_files[$profile]}"
                echo "  Added TF: $tf_name (from $meme_name) to profile: $profile"
            else
                echo "  WARNING: No profile mapping found for $tf_name"
            fi
        else
            echo "  WARNING: Missing file for $meme_name"
        fi
    done
    
    # Create output directory
    mkdir -p "${PLOTS_DIR}"
    
    # Create and submit jobs for each profile
    for profile in "${PROFILES[@]}"; do
        tf_list_file="${profile_files[$profile]}"
        tf_count=$(wc -l < "$tf_list_file")
        
        echo ""
        echo "Profile: $profile - Total TFs: $tf_count"
        
        # Skip if no TFs found for this profile
        if [[ $tf_count -eq 0 ]]; then
            echo "No TFs found for $cell_line - $profile, skipping..."
            continue
        fi
        
        # Generate the SLURM job script for this cell line and profile
        job_script="${JOB_SCRIPTS_DIR}/graphdnpsnucocc_${cell_line}_${profile}.sh"
        
        cat > "$job_script" <<'JOBSCRIPT'
#!/bin/bash -l
#SBATCH -J CELL_PROFILE
#SBATCH -t 06:00:00
#SBATCH --mail-user=PLACEHOLDER_EMAIL
#SBATCH -A PLACEHOLDER_ACCOUNT -p PLACEHOLDER_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40g
#SBATCH -o OUTPUT_FILE
#SBATCH -e ERROR_FILE
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Run the averaging script
python SCRIPTS_DIR/makeavgdnpsprofile.py \
    "TF_LIST_FILE" \
    "OUTPUT_PNG" \
    "CELL_LINE" \
    "PROFILE"

echo "Processed TF_COUNT TFs for CELL_LINE - PROFILE"
JOBSCRIPT
        
        # Replace placeholders
        sed -i "s|CELL_PROFILE|${cell_line}_${profile}|g" "$job_script"
        sed -i "s|OUTPUT_FILE|${JOB_SCRIPTS_DIR}/${cell_line}_${profile}_output.txt|g" "$job_script"
        sed -i "s|ERROR_FILE|${JOB_SCRIPTS_DIR}/${cell_line}_${profile}_error.txt|g" "$job_script"
        sed -i "s|TF_LIST_FILE|${tf_list_file}|g" "$job_script"
        sed -i "s|OUTPUT_PNG|${PLOTS_DIR}/${cell_line}_hg38_${profile}_average_nsm_occ_deltaNPS_profile.png|g" "$job_script"
        sed -i "s|CELL_LINE|${cell_line}|g" "$job_script"
        sed -i "s|PROFILE|${profile}|g" "$job_script"
        sed -i "s|TF_COUNT|${tf_count}|g" "$job_script"
        sed -i "s|SCRIPTS_DIR|${SCRIPTS_DIR}|g" "$job_script"
        sed -i "s|PLACEHOLDER_EMAIL|${EMAIL}|g" "$job_script"
        sed -i "s|PLACEHOLDER_ACCOUNT|${ACCOUNT}|g" "$job_script"
        sed -i "s|PLACEHOLDER_PARTITION|${PARTITION}|g" "$job_script"
        
        # Make the job script executable
        chmod +x "$job_script"
        
        # Submit the job
        echo "Submitting job for $cell_line - $profile with $tf_count TFs"
        sbatch "$job_script"
    done
    
    echo ""
    echo "================================"
done

echo "All jobs submitted!"
