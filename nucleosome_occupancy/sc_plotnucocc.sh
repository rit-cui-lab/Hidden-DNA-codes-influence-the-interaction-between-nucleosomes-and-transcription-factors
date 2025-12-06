#!/bin/bash
# sc_plotnucocc.sh
# Author: Chris Carson
# Description: Generates SLURM jobs for plotting nucleosome occupancy profiles.
# Scans directories for profile files and creates plotting jobs.

BASE_DIR="./chipseq/nucleosome_occupancy"
OUTPUT_DIR="./plot_nuc_occ_slurm_scripts"
PYTHON_SCRIPT="./scripts/nucleosome_occupancy/plotnucocc.py"

CELL_LINES=("H1_hESC")

mkdir -p "$OUTPUT_DIR"

job_count=0

for cell_line in "${CELL_LINES[@]}"; do
    cell_line_dir="$BASE_DIR/$cell_line"
    
    if [ ! -d "$cell_line_dir" ]; then
        echo "Warning: Cell line directory not found: $cell_line_dir"
        continue
    fi
    
    echo "Processing cell line: $cell_line"
    
    for meme_dir in "$cell_line_dir"/*; do
        [ ! -d "$meme_dir" ] && continue
        
        meme_name=$(basename "$meme_dir")
        profile_file="$meme_dir/${meme_name}_best_score_nsm_profile.tab"
        
        if [ -f "$profile_file" ]; then
            slurm_script="$OUTPUT_DIR/${cell_line}_${meme_name}_job.sh"
            
            cat > "$slurm_script" << EOF
#!/bin/bash -l
#SBATCH -t 00:30:00
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=4g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --job-name=${cell_line}_${meme_name}
#SBATCH --output=${cell_line}_${meme_name}_%j.out
#SBATCH --error=${cell_line}_${meme_name}_%j.err

python ${PYTHON_SCRIPT} -d "$profile_file" -c "$cell_line" -o ./results/
EOF
            
            chmod +x "$slurm_script"
            echo "Created Slurm script: $slurm_script"
            ((job_count++))
        fi
    done
done

echo "========================================="
echo "Total job scripts created: $job_count"
echo "Scripts located in: $OUTPUT_DIR"
echo "========================================="

# Create submission script
submit_script="$OUTPUT_DIR/submit_all_jobs.sh"
cat > "$submit_script" << 'EOF'
#!/bin/bash
for script in "$OUTPUT_DIR"/*_job.sh; do
    echo "Submitting $script"
    sbatch "$script"
done
EOF

chmod +x "$submit_script"
echo "Created submission script: $submit_script"
