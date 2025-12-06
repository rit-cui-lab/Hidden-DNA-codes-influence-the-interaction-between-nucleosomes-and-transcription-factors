#!/usr/bin/env python3
# calculate_GC_content.py
# Author: Chris Carson, Feng Cui
# Description: Calculates GC content percentage from a FASTA file.
# Counts G and C nucleotides and outputs final GC percentage.

import sys

def GC_content(seq):
    """Count G and C nucleotides in sequence"""
    GC_count = 0
    for char in seq:
        if char == 'G' or char == 'C':
            GC_count += 1
    return GC_count

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    GC = 0
    total_length = 0
    
    with open(input_file, 'r') as ih, open(output_file, 'w') as oh:
        for line in ih.readlines():
            if line.startswith('>'):
                name = line.strip()
            else:
                seq = line.strip().upper()
                GC += GC_content(seq)
                total_length += len(seq)
        
        final_GC = GC / total_length
        entry = "{:.2f}\n".format(final_GC)
        oh.write(entry)
