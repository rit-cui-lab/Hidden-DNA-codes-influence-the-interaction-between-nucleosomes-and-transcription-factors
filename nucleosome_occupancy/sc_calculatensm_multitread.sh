#!/bin/bash
# sc_calculatensm_multithread.sh
# Author: Chris Carson
# Description: Submits batch SLURM jobs for multithread NSM calculation.
# Processes chromosome files in parallel with configurable resources.

MAX_CORES=18
MEMORY="128G"
TIME="12:00:00"

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/directory"
    echo "Example: $0 ./mnase/chromosome_files/"
    exit 1
fi

if [ ! -d "$1" ]; then
    echo "Error: $1 is not a directory"
    exit 1
fi

# Process each chromosome file
for chr_file in "$1"/*chr*; do
    if [ ! -f "$chr_file" ]; then
        continue
    fi
    
    # Extract chromosome number from filename
    chr_name=$(basename "$chr_file" | sed 's/.*chr\([^_]*\).*/\1/')
    
    # Create SLURM submission
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=nuc_occ_chr${chr_name}
#SBATCH --output=nsm_chr${chr_name}_%j.out
#SBATCH --error=nsm_chr${chr_name}_%j.err
#SBATCH --cpus-per-task=${MAX_CORES}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${TIME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -A factor -p tier3
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

python ./scripts/nucleosome_occupancy/calculatensm_multithread.py ${chr_file}
EOF
    
    echo "Submitted job for chromosome ${chr_name}"
done
