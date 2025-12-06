#!/bin/bash
# =============================================================================
# callpipeline2.sh
# Author: Chris Carson
#
# Summary:
# This script generates and submits SLURM array jobs for Part 2 of the dNPS pipeline.
# It intersects the binned TF binding site regions with MNase-seq nucleosome dyad
# positions to identify nucleosomes near TF binding sites. Uses a multi-node array
# job configuration for parallel processing across 200 iterations (bins).
#
# The script:
#   1. Creates worker scripts that process subsets of iterations
#   2. Uses SLURM arrays to distribute work across multiple nodes
#   3. Intersects binned BED files with chromosome-split MNase bedgraph files
#   4. Combines results into final output files
#
# Prerequisites:
#   - Part 1 output (binned BED files)
#   - MNase-seq dyad bedgraph files split by chromosome
#   - bedtools installed
#
# Output:
#   - final_best_site_sorted_unique_{1-200}.bed files per TF
# =============================================================================

# =============================================================================
# CONFIGURATION - Modify these paths for your environment
# =============================================================================
BASE_DIR="/path/to/your/project"
DNPS_DIR="${BASE_DIR}/chipseq/dNPS"
NUCOCC_DIR="${BASE_DIR}/chipseq/nucleosome_occupancy"
MNASE_DIR="${BASE_DIR}/mnase"
SCRIPTS_DIR="${BASE_DIR}/scripts/dnps"
JOB_SCRIPTS_DIR="${SCRIPTS_DIR}/pt2_job_scripts"
EMAIL="your_email@example.com"
ACCOUNT="your_account"
PARTITION="your_partition"
SPACK_ENV="your_spack_env"
# =============================================================================

mkdir -p "$JOB_SCRIPTS_DIR"

CELL_LINES=("H1_hESC")

CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

# Multi-node configuration
NODES_PER_JOB=2
ITERATIONS_PER_ARRAY=20
TOTAL_ITERATIONS=200
NUM_ARRAY_JOBS=$((TOTAL_ITERATIONS / ITERATIONS_PER_ARRAY))

echo "Configuration:"
echo "  - Total iterations: $TOTAL_ITERATIONS"
echo "  - Array jobs: $NUM_ARRAY_JOBS"
echo "  - Iterations per array job: $ITERATIONS_PER_ARRAY"
echo "  - Nodes per array job: $NODES_PER_JOB"
echo "  - Total node-hours: $((NUM_ARRAY_JOBS * NODES_PER_JOB))"

for cell_line in "${CELL_LINES[@]}"; do
    echo "Processing cell line: $cell_line"

    for meme_dir in "$DNPS_DIR/$cell_line"/*.meme/; do
        [[ ! -d "$meme_dir" ]] && break

        meme_name=$(basename "$meme_dir")
        echo "  Processing meme directory: $meme_name"

        bed_file="${meme_dir}best_site_sorted_unique.bed"

        echo "    Creating multi-node array job for: $bed_file"

        output_dir="$meme_dir"
        
        echo "    Checking for input files..."
        sample_file="${output_dir}/${meme_name}_${cell_line}_best_site_sorted_unique_1.bed"
        if [[ ! -f "$sample_file" ]]; then
            echo "    ERROR: Input files not found. Expected pattern: ${meme_name}_${cell_line}_best_site_sorted_unique_*.bed"
            continue
        fi
        echo "    Input files found - proceeding with job creation"

        unique_id="${cell_line}_${meme_name}_multinode_array"
        
        worker_script="${JOB_SCRIPTS_DIR}/worker_${unique_id}.sh"
        
cat > "$worker_script" << 'EOF'
#!/bin/bash

# Get task parameters
ARRAY_TASK_ID=$1
TASK_ID=$SLURM_PROCID
HOSTNAME=$(hostname)

echo "=== Task Info ==="
echo "Array Task: $ARRAY_TASK_ID"
echo "Task ID: $TASK_ID"
echo "Hostname: $HOSTNAME"
echo "================="

# Hard-coded values
ITERATIONS_PER_ARRAY=20
NODES_PER_JOB=2

# Calculate iteration range for this array job
START_ITER=$(( ($ARRAY_TASK_ID - 1) * $ITERATIONS_PER_ARRAY + 1 ))
END_ITER=$(( $ARRAY_TASK_ID * $ITERATIONS_PER_ARRAY ))

# Calculate iterations per task
ITERATIONS_PER_TASK=$(( $ITERATIONS_PER_ARRAY / $NODES_PER_JOB ))

# Calculate iteration ranges for this task
if [ $TASK_ID -eq 0 ]; then
    TASK_START=$START_ITER
    TASK_END=$(( $START_ITER + $ITERATIONS_PER_TASK - 1 ))
else
    TASK_START=$(( $START_ITER + $ITERATIONS_PER_TASK ))
    TASK_END=$END_ITER
fi

echo "Task $TASK_ID on $HOSTNAME processing iterations $TASK_START to $TASK_END"

# Change to work directory
cd "PLACEHOLDER_OUTPUT_DIR"

# Process assigned iterations
for i in $(seq $TASK_START $TASK_END); do
    echo "Task $TASK_ID: Processing iteration $i at $(date)"
    
    # Process each chromosome
PLACEHOLDER_CHROMOSOME_COMMANDS
    
    # Combine results
    cat chr*_iter${i}.bed > result_${i}.bed
    mv result_${i}.bed final_best_site_sorted_unique_${i}.bed
    rm -f chr*_iter${i}.bed
    
    echo "Task $TASK_ID: Completed iteration $i"
done

echo "Task $TASK_ID completed all iterations"
EOF

        # Replace placeholders in worker script
        sed -i "s|PLACEHOLDER_OUTPUT_DIR|${output_dir}|g" "$worker_script"
        
        # Build chromosome commands
        chr_commands=""
        for chr in "${CHROMOSOMES[@]}"; do
            chr_commands="${chr_commands}    bedtools intersect -a ${output_dir}/${meme_name}_${cell_line}_best_site_sorted_unique_\${i}.bed -b ${MNASE_DIR}/${cell_line}/bedgraphs/split_filtered_${cell_line}_dyads_final_chr${chr}.bg > chr${chr}_iter\${i}.bed"$'\n'
        done
        
        # Replace chromosome commands placeholder
        temp_file=$(mktemp)
        echo "$chr_commands" > "$temp_file"
        sed -i "/PLACEHOLDER_CHROMOSOME_COMMANDS/r $temp_file" "$worker_script"
        sed -i "/PLACEHOLDER_CHROMOSOME_COMMANDS/d" "$worker_script"
        rm "$temp_file"

        chmod +x "$worker_script"

        # Create the main sbatch script
        job_script="${JOB_SCRIPTS_DIR}/part2_${unique_id}.sh"
        
cat > "$job_script" << 'EOF'
#!/bin/bash
#SBATCH --job-name=PLACEHOLDER_JOB_NAME_%A_%a
#SBATCH --time=05:00:00
#SBATCH --mail-user=PLACEHOLDER_EMAIL
#SBATCH --account=PLACEHOLDER_ACCOUNT
#SBATCH --partition=PLACEHOLDER_PARTITION
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=35g
#SBATCH --array=1-10
#SBATCH --output=PLACEHOLDER_OUTPUT_PREFIX_%A_%a_output.txt
#SBATCH --error=PLACEHOLDER_OUTPUT_PREFIX_%A_%a_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Load environment
spack env activate PLACEHOLDER_SPACK_ENV

echo "=== SLURM Job Info ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "Number of tasks: $SLURM_NTASKS"
echo "Nodelist: $SLURM_JOB_NODELIST"
echo "======================"

# Use srun to launch the worker script on all allocated tasks
srun PLACEHOLDER_WORKER_SCRIPT $SLURM_ARRAY_TASK_ID

echo "Array job $SLURM_ARRAY_TASK_ID completed"
EOF

        # Replace placeholders in job script
        sed -i "s|PLACEHOLDER_JOB_NAME|${unique_id}|g" "$job_script"
        sed -i "s|PLACEHOLDER_OUTPUT_PREFIX|${JOB_SCRIPTS_DIR}/${unique_id}|g" "$job_script"
        sed -i "s|PLACEHOLDER_WORKER_SCRIPT|${worker_script}|g" "$job_script"
        sed -i "s|PLACEHOLDER_EMAIL|${EMAIL}|g" "$job_script"
        sed -i "s|PLACEHOLDER_ACCOUNT|${ACCOUNT}|g" "$job_script"
        sed -i "s|PLACEHOLDER_PARTITION|${PARTITION}|g" "$job_script"
        sed -i "s|PLACEHOLDER_SPACK_ENV|${SPACK_ENV}|g" "$job_script"
        
        chmod +x "$job_script"
        
        echo "    Created sbatch script: $job_script"
        echo "    Created worker script: $worker_script"
        
    done
    echo "Completed processing cell line: $cell_line"
done
