#!/usr/bin/env python3
# plotnucocc.py
# Author: Chris Carson
# Description: Creates publication-quality nucleosome occupancy profile plots.
# Parses data files, applies smoothing, generates symmetrical profiles.

import matplotlib.pyplot as plt
import numpy as np
import re
import sys
import argparse
import os

def parse_filename(filename):
    """Extract transcription factor and database info from filename"""
    basename = os.path.basename(filename)
    parts = basename.split('_')
    transcription_factor = parts[0]
    
    memefile_db = None
    for part in parts[1:]:
        if '.meme' in part:
            memefile_db = part.split('.meme')[0]
            break
    
    return transcription_factor, memefile_db

def read_data_file(filename, cell_line_name):
    """Read data from file, searching for provided cell line"""
    with open(filename, 'r') as file:
        content = file.read()
    
    # Extract bin labels to determine min and max
    match = re.search(r'bin labels\s+([-\d.]+\w+)\s+center\s+([-\d.]+\w+)', content)
    if match:
        min_label = match.group(1)
        max_label = match.group(2)
        min_val = float(min_label.replace('Kb', '')) * 1000 
        max_val = float(max_label.replace('Kb', '')) * 1000
    else:
        min_val = -1000
        max_val = 1000
    
    lines = content.strip().split('\n')
    data_lines = [line for line in lines if line.strip() and not line.startswith('bin labels')]
    
    # Find line with specified cell line
    cell_line_data_line = None
    for line in data_lines:
        if line.startswith(cell_line_name) or cell_line_name in line.split():
            cell_line_data_line = line
            break
    
    if cell_line_data_line is None:
        raise ValueError(f"Could not find cell line '{cell_line_name}' in file.")
    
    # Extract values
    data_parts = cell_line_data_line.split()
    cell_line_values = [float(val) for val in data_parts[1:] if is_number(val)]
    
    # Create x-axis values based on range and number of points
    num_points = len(cell_line_values)
    x_values = np.linspace(min_val, max_val, num_points)
    
    return x_values, cell_line_values

def smooth_data(data_list, window_size=7):
    """Apply sliding window smoothing"""
    if len(data_list) < window_size:
        return data_list.copy()
    
    smoothed_data = []
    
    # Beginning: use mean of first window
    first_window_mean = np.mean(data_list[0:window_size])
    for i in range(window_size // 2):
        smoothed_data.append(first_window_mean)
    
    # Middle: sliding window
    for i in range(len(data_list) - window_size + 1):
        avg = np.mean(data_list[i:i + window_size])
        smoothed_data.append(avg)
    
    # End: use mean of last window
    last_window_mean = np.mean(data_list[-window_size:])
    for i in range(window_size // 2):
        smoothed_data.append(last_window_mean)
    
    return smoothed_data

def is_number(s):
    """Check if string can be converted to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False

def make_symmetrical(x_values, y_values):
    """Create symmetrical profile by averaging equidistant values"""
    center_idx = len(x_values) // 2
    sym_y_values = np.zeros_like(y_values)
    
    for i in range(len(y_values)):
        mirror_idx = 2 * center_idx - i
        
        if 0 <= mirror_idx < len(y_values):
            sym_y_values[i] = (y_values[i] + y_values[mirror_idx]) / 2
        else:
            sym_y_values[i] = y_values[i] / 2
    
    return sym_y_values

def plot_data(filename, cell_line_name, symmetrical=True, smooth=True, output_dir="./"):
    """Plot nucleosome occupancy profile"""
    try:
        tf, memefile = parse_filename(filename)
        
        if tf and memefile:
            title_prefix = f"{tf}_{memefile}"
        else:
            title_prefix = os.path.basename(filename).split('.')[0]
        
        x_values, y_values = read_data_file(filename, cell_line_name)
        
        # Remove rightmost point and recalculate x range
        y_values = y_values[:-1]
        num_points = len(y_values)
        x_min = x_values[0]
        x_max = x_values[-1]
        x_values = np.linspace(x_min, x_max, num_points)
        
        # Apply smoothing if requested
        if smooth:
            y_values = smooth_data(y_values)
        
        # Set y-axis limits
        y_min = min(y_values) - 0.05
        y_max = max(y_values) + 0.05
        
        # X-axis ticks every 100 bp
        tick_positions = np.arange(-1000, 1001, 100)
        
        plt.figure(figsize=(10, 6))
        
        if symmetrical:
            y_values_sym = make_symmetrical(x_values, y_values)
            sym_y_min = min(y_values_sym) - 0.05
            sym_y_max = max(y_values_sym) + 0.05
            
            plt.plot(x_values, y_values_sym, 'r-', linewidth=1.5)
            plt.title(f'{title_prefix} - Symmetrical {cell_line_name} Nucleosome Occupancy')
            plt.ylim(sym_y_min, sym_y_max)
        else:
            plt.plot(x_values, y_values, 'k-', linewidth=1.5)
            plt.title(f'{title_prefix} - {cell_line_name} Nucleosome Occupancy')
            plt.ylim(y_min, y_max)
        
        plt.xlabel('Distance from center (bp)')
        plt.ylabel('Nucleosome Occupancy')
        plt.grid(True, alpha=0.3)
        plt.xticks(tick_positions, [str(x) for x in tick_positions])
        plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        plt.tight_layout()
        
        # Create output filename
        suffix_parts = []
        if symmetrical:
            suffix_parts.append("symmetrical")
        else:
            suffix_parts.append("original")
        if smooth:
            suffix_parts.append("smoothed")
        
        suffix = "_".join(suffix_parts)
        output_filename = f"{title_prefix}_{cell_line_name}_{suffix}.png"
        output_path = os.path.join(output_dir, output_filename)
        
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_path}")
        plt.close()
        
        return True
        
    except Exception as e:
        print(f"Error processing {filename} for {cell_line_name}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Plot nucleosome occupancy profile.')
    parser.add_argument('-d', dest='filename', required=True, help='Input data file')
    parser.add_argument('-c', dest='cell_line', required=True, help='Cell line name to plot')
    parser.add_argument('-o', dest='output_dir', default="./", help='Output directory')
    parser.add_argument('--original', dest='symmetrical', action='store_false', 
                        help='Plot original data instead of symmetrical')
    parser.add_argument('--no-smooth', dest='smooth', action='store_false',
                        help='Disable smoothing')
    parser.set_defaults(symmetrical=True, smooth=True)
    
    args = parser.parse_args()
    
    result = plot_data(args.filename, args.cell_line, args.symmetrical, args.smooth, args.output_dir)
    if not result:
        print(f"Could not plot data for cell line '{args.cell_line}'")

if __name__ == "__main__":
    main()
