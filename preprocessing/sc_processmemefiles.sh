#!/bin/bash
# slurmcall_processmemefiles_all.sh
# Author: Chris Carson
# Description: Batch job submission script for processing multiple MEME files.
# Generates individual SLURM jobs for each MEME file in a directory.

if [ $# -ne 3 ]; then
    echo "Usage: $0 <meme_files_directory> <protein_file_path> <cell_line>"
    exit 1
fi

MEME_DIR="$1"
PROTEIN_PATH="$2"
CELL_LINE="$3"
OUTPUT_DIR="meme_job_scripts"
RESULTS_DIR="meme_results"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$RESULTS_DIR"

for meme_file in "$MEME_DIR"/*.meme; do
    [ -e "$meme_file" ] || continue
    
    base_name=$(basename "$meme_file" .meme)
    
    mkdir -p "$RESULTS_DIR/$base_name"
    
    cat > "$OUTPUT_DIR/${base_name}_job.sh" << EOL
#!/bin/bash -l
#SBATCH -J ${base_name}
#SBATCH -t 2-7:0:0
#SBATCH -A factor -p tier3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100g
#SBATCH -o $RESULTS_DIR/${base_name}/${base_name}_output.txt
#SBATCH -e $RESULTS_DIR/${base_name}/${base_name}_error.txt
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd ./chipseq/memefiles/$CELL_LINE

python ./scripts/preprocessing/processmemefiles.py "$meme_file" "$PROTEIN_PATH"
EOL

    chmod +x "$OUTPUT_DIR/${base_name}_job.sh"
    
    echo "Created job script for $base_name"
done

echo "Done! Job scripts are in $OUTPUT_DIR"
echo "Results will be stored in $RESULTS_DIR"
