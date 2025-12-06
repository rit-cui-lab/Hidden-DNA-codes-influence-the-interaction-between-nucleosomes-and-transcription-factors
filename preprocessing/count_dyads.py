#!/usr/bin/env python3
#Author: ChrisCarson

import sys

def numeric_sort(key):
    return int(key)

# Generate output filename by adding _final before .txt
input_filename = sys.argv[1]
output_filename = input_filename.replace('.txt', '_final.txt')

with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output1:
    freq_start = {}
    
    # Process input file
    for line in input_file:
        line = line.strip()
        data = line.split('\t')     
        chr_val = data[0]
        dyad = data[1]
        
        # Create a tuple key to store both chromosome and dyad
        key = (chr_val, dyad)
        freq_start[key] = freq_start.get(key, 0) + 1    
    
    # Write output, unpacking chromosome and dyad from tuple key
    for (chr_val, dyad), count in sorted(freq_start.items(), key=lambda x: numeric_sort(x[0][1])):
        output1.write(f"{chr_val}\t{dyad}\t{count}\n")
