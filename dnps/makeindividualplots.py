#!/usr/bin/env python3
"""
makeindividualplots.py
Author: Chris Carson

Summary:
This script generates individual nucleosome positioning profile plots for single
transcription factors. Each plot shows three data series: in vivo nucleosome
occupancy (red), delta-NPS percentage (blue), and optionally in vitro nucleosome
occupancy (green) if available. The script handles missing in vitro data gracefully,
creating plots without the green line and logging skipped TFs to a tracking file.

Usage:
    python makeindividualplots.py <invivo_file> <meme_name> <output.png> <cell_line> <profile_type> <skipped_file>

Arguments:
    invivo_file   - Path to in vivo nucleosome occupancy file
    meme_name     - Full meme name (e.g., JUN_HOCOMOCOv11)
    output.png    - Output plot filename
    cell_line     - Cell line identifier
    profile_type  - Profile category (DIP_DIP, DIP_PEAK, PEAK_DIP, PEAK, AMBIGUOUS)
    skipped_file  - File to record TFs without in vitro data
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
    data_list = data_array.tolist() if isinstance(data_array, np.ndarray) else data_array
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


def make_plot(arr, arr2, arr3, output, cell_line, tf_name, profile_type, has_invitro):
    """Generate three-axis plot for individual TF"""
    lst = np.array(list(map(float, arr)))
    rev_lst = np.array(list(reversed(lst)))
    final_lst = 0.5 * lst + 0.5 * rev_lst
    
    lst_2 = np.array(arr2)
    rev_lst_2 = np.array(list(reversed(lst_2)))
    final_lst_2 = 0.5 * lst_2 + 0.5 * rev_lst_2
    
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
    
    # Third y-axis (right side 2) - In vitro nucleosome occupancy (green) - only if has data
    if has_invitro:
        lst_3 = np.array(list(map(float, arr3)))
        rev_lst_3 = np.array(list(reversed(lst_3)))
        final_lst_3 = 0.5 * lst_3 + 0.5 * rev_lst_3
        
        ax3 = ax1.twinx()
        ax3.spines['right'].set_position(('outward', 60))
        ax3.set_ylabel('In vitro nucleosome occupancy', color='green', fontsize=18)
        ax3.plot(x, final_lst_3, linewidth=2.5, color='green', label='In vitro')
        ax3.tick_params(axis='y', labelsize=16, labelcolor='green')
    
    # Convert DIP_DIP back to DIP/DIP for display
    profile_display = profile_type.replace('_', '/')
    
    # Extract clean TF name from meme_name
    clean_tf_name = tf_name.replace('_HOCOMOCOv11', '').replace('_complete-factorbook-catalog', '')
    
    # Update title
    title = f'{cell_line} - {clean_tf_name} ({profile_display})'
    if not has_invitro:
        title += ' [no invitro]'
    
    plt.title(title, fontsize=20)
    plt.margins(x=0)
    plt.tight_layout(pad=2.0)
    plt.savefig(output)
    plt.close()


if __name__ == "__main__":
    invivo_file = sys.argv[1]   # Path to in vivo nucleosome occupancy file
    meme_name = sys.argv[2]     # Full meme name (e.g., JUN_HOCOMOCOv11)
    output = sys.argv[3]        # Output plot path
    cell = sys.argv[4]          # Cell line
    profile_type = sys.argv[5]  # Profile type
    skipped_file = sys.argv[6]  # File to record skipped TFs
    
    # Read in vivo nucleosome occupancy
    with open(invivo_file, 'r') as ih:
        for line in ih.readlines()[2:]:
            data = line.strip().split()
            nsm_occ = list(map(float, data[1:]))
    
    # Read in vitro nucleosome occupancy
    dir_path = os.path.dirname(invivo_file)
    filename = os.path.basename(invivo_file)
    invitro_filename = "invitro_" + filename
    invitro_file = os.path.join(dir_path, invitro_filename)
    
    has_invitro = False
    invitro_occ = None
    
    # Check if invitro file exists and has data
    if not os.path.exists(invitro_file):
        clean_tf_name = meme_name.replace('_HOCOMOCOv11', '').replace('_complete-factorbook-catalog', '')
        with open(skipped_file, 'a') as f:
            f.write(f"{clean_tf_name}\t(file not found)\n")
        print(f"No invitro data for {meme_name} - creating plot without invitro line")
    elif os.path.getsize(invitro_file) == 0:
        clean_tf_name = meme_name.replace('_HOCOMOCOv11', '').replace('_complete-factorbook-catalog', '')
        with open(skipped_file, 'a') as f:
            f.write(f"{clean_tf_name}\t(empty file)\n")
        print(f"Empty invitro file for {meme_name} - creating plot without invitro line")
    else:
        # Try to read invitro data
        try:
            with open(invitro_file, 'r') as ih:
                lines = ih.readlines()
                # Check if there's actual data (more than just headers)
                if len(lines) <= 2:
                    clean_tf_name = meme_name.replace('_HOCOMOCOv11', '').replace('_complete-factorbook-catalog', '')
                    with open(skipped_file, 'a') as f:
                        f.write(f"{clean_tf_name}\t(no data rows)\n")
                    print(f"No data in invitro file for {meme_name} - creating plot without invitro line")
                else:
                    # Parse the data
                    for line in lines[2:]:
                        data = line.strip().split()
                        if data and len(data) >= 2:
                            invitro_occ = list(map(float, data[1:]))
                            has_invitro = True
                            break
                    
                    if not has_invitro:
                        clean_tf_name = meme_name.replace('_HOCOMOCOv11', '').replace('_complete-factorbook-catalog', '')
                        with open(skipped_file, 'a') as f:
                            f.write(f"{clean_tf_name}\t(invalid data format)\n")
                        print(f"Invalid invitro data for {meme_name} - creating plot without invitro line")
        except Exception as e:
            clean_tf_name = meme_name.replace('_HOCOMOCOv11', '').replace('_complete-factorbook-catalog', '')
            with open(skipped_file, 'a') as f:
                f.write(f"{clean_tf_name}\t(read error: {str(e)})\n")
            print(f"Error reading invitro file for {meme_name} - creating plot without invitro line")
    
    # Calculate deltaNPS
    deltaNPS_lst = []
    for i in range(200):
        filename = f"{DNPS_DIR}/{cell}/{meme_name}.meme/best_site_sorted_unique_{i+1}_147bp_four_types.txt"
        
        try:
            deltaNPS = calc_delta_NPS(filename)
            deltaNPS_lst.append(deltaNPS)
        except FileNotFoundError:
            print(f"Warning: No file found for {meme_name} at position {i+1}")
            deltaNPS_lst.append(0)
    
    # Apply smoothing
    nsm_occ_smoothed = smooth_data(nsm_occ)
    deltaNPS_smoothed = smooth_data(deltaNPS_lst)
    
    # Only smooth invitro if we have data
    if has_invitro:
        invitro_occ_smoothed = smooth_data(invitro_occ)
    else:
        invitro_occ_smoothed = None
    
    # Make plot
    make_plot(nsm_occ_smoothed, deltaNPS_smoothed, invitro_occ_smoothed, 
              output, cell, meme_name, profile_type, has_invitro)
    
    print(f"Created plot for {meme_name}")
