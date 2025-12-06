#!/usr/bin/env python3
# normalize_bedgraph.py
# Author: Chris Carson
# Description: Normalizes bedGraph files by dividing occupancy values by average.
# Maintains bedGraph format (chromosome, start, end, value).

import sys

def calculate_average(input_file):
    """Calculate average occupancy value from bedGraph file"""
    values = []
    with open(input_file, 'r') as ih:
        for line in ih.readlines():
            parts = list(line.strip().split())
            if len(parts) >= 4:
                try:
                    values.append(float(parts[3]))
                except ValueError:
                    continue
    
    if not values:
        raise ValueError(f"No valid numeric values found in {input_file}")
    
    return sum(values) / len(values)

def normalize_bedgraph(input_file, output_file, avg):
    """Normalize bedGraph values by dividing by average"""
    with open(input_file, 'r') as ih, open(output_file, 'w') as oh:
        for line in ih.readlines():
            parts = list(line.strip().split())
            if len(parts) >= 4:
                try:
                    # Divide the value (4th column) by average
                    parts[3] = str(float(parts[3]) / avg)
                    oh.write('\t'.join(parts))
                    oh.write('\n')
                except (ValueError, IndexError):
                    # Skip malformed lines
                    continue

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python normalize_bedgraph.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        # Calculate average
        avg = calculate_average(input_file)
        print(f"Average occupancy value: {avg:.6f}")
        
        # Normalize and write
        normalize_bedgraph(input_file, output_file, avg)
        print(f"Normalized bedGraph written to: {output_file}")
        
    except FileNotFoundError:
        print(f"Error: File {input_file} not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
