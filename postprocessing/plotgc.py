#!/usr/bin/env python3
# plotgc.py
# Author: Chris Carson
# Description: Creates publication-quality GC content plots.
# Applies sliding window smoothing and generates symmetric plots.

import sys
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean

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

def make_plot(data, output_file, protein_name):
    """Generate GC content plot with symmetry"""
    lst = np.array(list(map(float, data)))
    rev_lst = np.array(list(reversed(lst)))
    final_lst = 0.5 * lst + 0.5 * rev_lst
    x = np.array(range(len(final_lst)))
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(left=0.15, right=0.9, bottom=0.2)
    
    ax1.set_xlabel('Distance to TF center', color='black', fontsize=14)
    ax1.set_ylabel('GC content', color='black', fontsize=14)
    
    # Adjust x-axis labels based on data length
    if len(final_lst) == 200:
        plt.xticks([0, 19.5, 39.5, 59.5, 79.5, 99.5, 119.5, 139.5, 159.5, 179.5, 199], 
                   ['-1.0Kb', '', '', '', '', 'center', '', '', '', '', '1.0Kb'], fontsize=14)
    else:
        # Auto-scale for different lengths
        mid = len(final_lst) // 2
        plt.xticks([0, mid, len(final_lst)-1], 
                   ['-1.0Kb', 'center', '1.0Kb'], fontsize=14)
    
    ax1.tick_params(axis='y', which='major', labelsize=12)
    ax1.yaxis.set_ticks_position('both')
    
    line1, = ax1.plot(x, final_lst, linewidth=1.5, color='black', label="GC")
    
    plt.legend(frameon=False, fontsize=12, loc="upper right")
    plt.title(f"GC content - {protein_name}", fontsize=16)
    plt.margins(x=0)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python plotgc.py <input_file> <output_plot> <protein_name>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    protein_name = sys.argv[3]
    
    gc_values = []
    
    try:
        with open(input_file, 'r') as ih:
            for line in ih:
                line = line.strip()
                if line:  # Skip empty lines
                    data = line.split()
                    if data:
                        gc_values.append(float(data[0]))
        
        if not gc_values:
            print(f"Error: No data found in {input_file}")
            sys.exit(1)
        
        print(f"Read {len(gc_values)} GC values from {input_file}")
        
        smoothed_data = smooth_data(gc_values)
        make_plot(smoothed_data, output_file, protein_name)
        
    except FileNotFoundError:
        print(f"Error: File {input_file} not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
