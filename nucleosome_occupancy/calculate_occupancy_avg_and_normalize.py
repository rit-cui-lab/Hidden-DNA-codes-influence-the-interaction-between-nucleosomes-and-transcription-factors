#!/usr/bin/env python3
# calculate_occupancy_avg_and_normalize.py
# Author: Chris Carson
# Description: Normalizes nucleosome occupancy by average value.
# Calculates average occupancy and divides all scores by that average.

import sys

def calculate_average(input_file):
    """Calculate average occupancy value from file"""
    occ = []
    with open(input_file, 'r') as ih:
        for line in ih.readlines():
            lst = list(line.strip().split())
            occ.append(float(lst[3]))
    return sum(occ)/len(occ)

def process_file(input_file, output_file, avg):
    """Normalize occupancy values by average"""
    with open(input_file, 'r') as ih, open(output_file, 'w') as oh:
        for line in ih.readlines():
            lst = list(line.strip().split())
            lst[3] = str(float(lst[3])/avg)
            oh.write('\t'.join(lst))
            oh.write('\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # First calculate the average
    avg = calculate_average(input_file)
    print(f"Average value: {avg}")
    
    # Then process the file using the calculated average
    process_file(input_file, output_file, avg)
