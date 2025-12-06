#!/usr/bin/env python3
"""
makeavgdnpsprofile.py
Author: Chris Carson

Summary:
This script generates averaged nucleosome positioning profiles for groups of
transcription factors. It combines in vivo nucleosome occupancy, in vitro
nucleosome occupancy, and delta-NPS (nucleosome positioning score) data across
multiple TFs belonging to the same profile category (e.g., DIP/DIP, PEAK, etc.).
The output is a three-axis plot showing the relationship between intrinsic
nucleosome preferences and in vivo positioning around TF binding sites.

Usage:
    python makeavgdnpsprofile.py <tf_list_file> <output.png> <cell_line> <profile_type>

Arguments:
    tf_list_file  - File containing paths to TF data files and meme names
    output.png    - Output plot filename
    cell_line     - Cell line identifier
    profile_type  - Profile category (DIP_DIP, DIP_PEAK, PEAK_DIP, PEAK, AMBIGUOUS)
"""

import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean

# =============================================================================
# CONFIGURATION - Modify these paths for your environment
# =============================================================================
BASE_DIR = "/path/to/your/project"
DNPS_DIR = f"{BASE_DIR}/chipseq/dNPS"
# =============================================================================


def smooth_data(data_array):
    """Apply 7-element sliding window smoothing with edge padding"""
    data_list = data_array.tolist()
    smoothed = []
    
    # First 3 positions: use mean of first 7 elements
    for i in range(3):
        smoothed.append(mean(data_list[0:7]))
    
    # Middle positions: sliding window
    for i in range(194):
        avg_window = mean(data_list[i:i+7])
        smoothed.append(avg_window)
    
    # Last 3 positions: use mean of last 7 elements
    for i in range(3):
        smoothed.append(mean(data_list[193:200]))
    
    return np.array(smoothed)


def calc_delta_NPS(input):
    """Calculate delta-NPS from four nucleosome types"""
    with open(input, 'r') as ih:
        line = ih.readline()
        lst = list(line.strip().split('\t'))
        type1 = int(lst[0])
        type2 = int(lst[1])
        type3 = int(lst[2])
        type4 = int(lst[3])
        sum = type1 + type2 + type3 + type4
        return (type1 - type4) / sum * 100


def make_plot(arr, arr2, arr3, output, cell_line, profile_type, num_tfs, invitro_count):
    """Generate three-axis plot with in vivo, delta-NPS, and in vitro data"""
    lst = np.array(list(map(float, arr)))
    rev_lst = np.array(list(reversed(lst)))
    final_lst = 0.5 * lst + 0.5 * rev_lst
    
    lst_2 = np.array(arr2)
    rev_lst_2 = np.array(list(reversed(lst_2)))
    final_lst_2 = 0.5 * lst_2 + 0.5 * rev_lst_2
    
    lst_3 = np.array(list(map(float, arr3)))
    rev_lst_3 = np.array(list(reversed(lst_3)))
    final_lst_3 = 0.5 * lst_3 + 0.5 * rev_lst_3
    
    x = np.array(range(200))
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # First y-axis (left) - In vivo nucleosome occupancy (red)
    ax1.set_xlabel('Distance to TF center', color='black', fontsize=18)
    ax1.set_ylabel('In vivo nucleosome occupancy', color='red', fontsize=18)
    plt.xticks([0, 19.5, 39.5, 59.5, 79.5, 99.5, 119.5, 139.5, 159.5, 179.5, 199], 
               ['-1.0Kb', '', '', '', '', 'center', '', '', '', '', '1.0Kb'], fontsize=18)
    ax1.tick_params(axis='y', which='major', labelsize=16, labelcolor='red')
    ax1.plot(x, final_lst, linewidth=2.5, color='red', label='In vivo')
    
    # Second y-axis (right side 1) - deltaNPS (blue)
    ax2 = ax1.twinx()
    ax2.set_ylabel('$\\Delta$NPS, %', color='blue', fontsize=18)
    ax2.plot(x, final_lst_2, linewidth=1.5, color='blue', label='$\\Delta$NPS')
    ax2.tick_params(axis='y', labelsize=16, labelcolor='blue')
    
    # Third y-axis (right side 2) - In vitro nucleosome occupancy (green)
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))
    ax3.set_ylabel('In vitro nucleosome occupancy', color='green', fontsize=18)
    ax3.plot(x, final_lst_3, linewidth=2.5, color='green', label='In vitro')
    ax3.tick_params(axis='y', labelsize=16, labelcolor='green')
    
    # Convert DIP_DIP back to DIP/DIP for display
    profile_display = profile_type.replace('_', '/')
    
    # Update title to show both total TFs and invitro count
    if invitro_count == num_tfs:
        title = f'{cell_line} - {profile_display} (n={num_tfs})'
    else:
        title = f'{cell_line} - {profile_display} (n={num_tfs}, invitro={invitro_count})'
    
    plt.title(title, fontsize=20)
    plt.margins(x=0)
    plt.tight_layout(pad=2.0)
    plt.savefig(output)


if __name__ == "__main__":
    tf_list_file = sys.argv[1]  # File containing list of TF info
    output = sys.argv[2]
    cell = sys.argv[3]
    profile_type = sys.argv[4]  # Profile type (DIP_DIP, DIP_PEAK, PEAK, or AMBIGUOUS)
    
    # Read the list of TFs and their input files
    tf_data = []
    with open(tf_list_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                input_file = parts[0]
                meme_name = parts[1]
                tf_data.append((input_file, meme_name))
    
    all_nsm_occ = []
    all_invitro_occ = []
    all_deltaNPS = []
    tfs_with_invitro = 0
    tfs_without_invitro = []
    
    # Process each TF
    for input_file, meme_name in tf_data:
        # Read in vivo nucleosome occupancy
        with open(input_file, 'r') as ih:
            for line in ih.readlines()[2:]:
                data = line.strip().split()
                nsm_occ = list(map(float, data[1:]))
        all_nsm_occ.append(nsm_occ)
        
        # Read in vitro nucleosome occupancy
        dir_path = os.path.dirname(input_file)
        filename = os.path.basename(input_file)
        invitro_filename = "invitro_" + filename
        invitro_file = os.path.join(dir_path, invitro_filename)
        
        try:
            with open(invitro_file, 'r') as ih:
                for line in ih.readlines()[2:]:
                    data = line.strip().split()
                    invitro_occ = list(map(float, data[1:]))
            all_invitro_occ.append(invitro_occ)
            tfs_with_invitro += 1
        except FileNotFoundError:
            print(f"Warning: No invitro file found for {meme_name}")
            tfs_without_invitro.append(meme_name)
        
        # Calculate deltaNPS for this TF
        deltaNPS_lst = []
        for i in range(200):
            filename = f"{DNPS_DIR}/{cell}/{meme_name}.meme/best_site_sorted_unique_{i+1}_147bp_four_types.txt"
            
            try:
                deltaNPS = calc_delta_NPS(filename)
                deltaNPS_lst.append(deltaNPS)
            except FileNotFoundError:
                print(f"Warning: No file found for {meme_name} at position {i+1}")
                deltaNPS_lst.append(0)
        
        all_deltaNPS.append(deltaNPS_lst)
    
    # Print summary of invitro data availability
    print(f"\n=== In vitro data summary ===")
    print(f"Total TFs: {len(tf_data)}")
    print(f"TFs with invitro data: {tfs_with_invitro}")
    print(f"TFs without invitro data: {len(tfs_without_invitro)}")
    if tfs_without_invitro:
        print(f"TFs missing invitro files: {', '.join(tfs_without_invitro)}")
    print(f"============================\n")
    
    # Calculate averages across all TFs
    avg_nsm_occ = np.mean(all_nsm_occ, axis=0)
    avg_deltaNPS = np.mean(all_deltaNPS, axis=0)
    
    # Calculate invitro average only if we have data
    if all_invitro_occ:
        avg_invitro_occ = np.mean(all_invitro_occ, axis=0)
    else:
        print("WARNING: No invitro data available for any TFs. Using zeros for invitro plot.")
        avg_invitro_occ = np.zeros(200)
    
    # Apply smoothing to all three datasets
    avg_nsm_occ_smoothed = smooth_data(avg_nsm_occ)
    avg_invitro_occ_smoothed = smooth_data(avg_invitro_occ)
    avg_deltaNPS_smoothed = smooth_data(avg_deltaNPS)
    
    # Save average deltaNPS (smoothed)
    outfile = f"{cell}_{profile_type}_average_deltaNPS_profile_smoothed.txt"
    with open(outfile, 'w') as oh:
        oh.write('\n'.join(str(i) for i in avg_deltaNPS_smoothed))
    
    # Make plot with smoothed in vivo, deltaNPS, and in vitro
    invitro_count = len(all_invitro_occ)
    make_plot(avg_nsm_occ_smoothed, avg_deltaNPS_smoothed, avg_invitro_occ_smoothed, 
              output, cell, profile_type, len(tf_data), invitro_count)
