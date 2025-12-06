#!/bin/bash
# =============================================================================
# callindividualplots.sh
# Author: Chris Carson
#
# Summary:
# This script generates and submits SLURM jobs to create individual nucleosome
# positioning profile plots for each transcription factor. Unlike the averaged
# plots, these show data for single TFs, allowing examination of TF-specific
# nucleosome positioning patterns.
#
# Each plot displays:
#   - In vivo nucleosome occupancy (red)
#   - Delta-NPS percentage (blue)
#   - In vitro nucleosome occupancy (green, if available)
#
# TFs are organized by profile category and cell line. The script handles
# missing in vitro data gracefully, creating plots without the green line
# and logging which TFs were affected.
#
# Prerequisites:
#   - Completed dNPS pipeline (parts 1-3)
#   - Nucleosome occupancy profile files
#   - TF-to-profile mapping files
#   - makeindividualplots.py script
#
# Output:
#   - Individual profile plots organized by cell line and profile category
# =============================================================================

# =============================================================================
# CONFIGURATION - Modify these paths for your environment
# =============================================================================
BASE_DIR="/path/to/your/project"
DNPS_DIR="${BASE_DIR}/chipseq/dNPS"
NUCOCC_DIR="${BASE_DIR}/chipseq/nucleosome_occupancy"
SCRIPTS_DIR="${BASE_DIR}/scripts/dnps"
JOB_SCRIPTS_DIR="${SCRIPTS_DIR}/plot_job_scripts/avgs_individual"
BASE_PLOTS_DIR="${DNPS_DIR}/plots/individual"
TF_PROFILE_MAPPING_DIR="${SCRIPTS_DIR}/avgs/tf_lists"
EMAIL="your_email@example.com"
ACCOUNT="your_account"
PARTITION="your_partition"
# =============================================================================

# Create directories if they don't exist
mkdir -p "$JOB_SCRIPTS_DIR"
mkdir -p "$BASE_PLOTS_DIR"

# Array of cell lines
CELL_LINES=("MCF-7" "H1_hESC" "HeLa" "IMR90" "HepG2" "HCT116")

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
    
    # Create output directories for each profile
    for profile in "${PROFILES[@]}"; do
        # Replace slashes with underscores for directory names
        profile_dir=$(echo "$profile" | tr '/' '_')
        profile_plots_dir="${BASE_PLOTS_DIR}/${cell_line}/${profile_dir}"
        mkdir -p "$profile_plots_dir"
    done
    
    # Collect all TFs for this cell line and categorize by profile
    declare -A profile_tfs
    for profile in "${PROFILES[@]}"; do
        profile_tfs[$profile]=""
    done
    
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
                # Add to this profile's list
                if [[ -z "${profile_tfs[$profile]}" ]]; then
                    profile_tfs[$profile]="$nucocc_file|$meme_name"
                else
                    profile_tfs[$profile]="${profile_tfs[$profile]} $nucocc_file|$meme_name"
                fi
                echo "  Added TF: $tf_name (from $meme_name) to profile: $profile"
            else
                echo "  WARNING: No profile mapping found for $tf_name"
            fi
        else
            echo "  WARNING: Missing file for $meme_name"
        fi
    done
    
    # Create and submit jobs for each profile
    for profile in "${PROFILES[@]}"; do
        tf_data="${profile_tfs[$profile]}"
        
        # Skip if no TFs found for this profile
        if [[ -z "$tf_data" ]]; then
            echo "  No TFs found for $cell_line - $profile, skipping..."
            continue
        fi
        
        # Count TFs
        tf_count=$(echo "$tf_data" | wc -w)
        
        echo ""
        echo "Profile: $profile - Creating individual plots for $tf_count TFs"
        
        # Replace slashes with underscores for directory names
        profile_dir=$(echo "$profile" | tr '/' '_')
        profile_plots_dir="${BASE_PLOTS_DIR}/${cell_line}/${profile_dir}"
        
        # Write TF data to a temporary file
        tf_data_file="${JOB_SCRIPTS_DIR}/${cell_line}_${profile}_tf_data.txt"
        echo "$tf_data" | tr ' ' '\n' > "$tf_data_file"
        
        # Create file to track TFs without invitro data
        skipped_file="${profile_plots_dir}/tfs_without_invitro.txt"
        > "$skipped_file"  # Initialize empty file
        
        # Generate the SLURM job script for this profile
        job_script="${JOB_SCRIPTS_DIR}/individual_plots_${cell_line}_${profile}.sh"
        
        cat > "$job_script" <<JOBSCRIPT
#!/bin/bash -l
#SBATCH -J ${cell_line}_${profile}_indiv
#SBATCH -t 06:00:00
#SBATCH --mail-user=${EMAIL}
#SBATCH -A ${ACCOUNT} -p ${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40g
#SBATCH -o ${JOB_SCRIPTS_DIR}/${cell_line}_${profile}_individual_output.txt
#SBATCH -e ${JOB_SCRIPTS_DIR}/${cell_line}_${profile}_individual_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Read TF data from file
while read -r tf_pair; do
    # Skip if empty
    [[ -z "\$tf_pair" ]] && continue
    
    # Split on pipe character
    invivo_file="\${tf_pair%%|*}"
    meme_name="\${tf_pair##*|}"
    
    # Create output filename
    clean_tf_name=\$(echo "\$meme_name" | sed 's/_HOCOMOCOv11$//' | sed 's/_complete-factorbook-catalog$//')
    output_file="${profile_plots_dir}/\${clean_tf_name}.png"
    
    echo "Processing \$clean_tf_name..."
    
    # Run the plotting script for this individual TF
    python ${SCRIPTS_DIR}/makeindividualplots.py \\
        "\$invivo_file" \\
        "\$meme_name" \\
        "\$output_file" \\
        "${cell_line}" \\
        "${profile}" \\
        "${skipped_file}"
        
done < "${tf_data_file}"

# Print summary of skipped TFs
if [[ -s "${skipped_file}" ]]; then
    skipped_count=\$(wc -l < "${skipped_file}")
    echo ""
    echo "=== Summary ==="
    echo "\$skipped_count TFs plotted without invitro data (missing or invalid files)"
    echo "See: ${skipped_file}"
else
    echo ""
    echo "=== Summary ==="
    echo "All TFs had valid invitro data"
    rm "${skipped_file}"
fi

echo "Completed all individual plots for ${cell_line} - ${profile}"
JOBSCRIPT
        
        # Make the job script executable
        chmod +x "$job_script"
        
        # Submit the job
        echo "  Submitting job for $cell_line - $profile with $tf_count TFs"
        sbatch "$job_script"
    done
    
    echo ""
    echo "================================"
done

echo "All individual plot jobs submitted!"
echo "Plots will be organized in: $BASE_PLOTS_DIR/{cell_line}/{profile}/"
