#!/usr/bin/env python3
# generate_10bp_bed.py
# Author: Chris Carson, Feng Cui
# Description: Converts dyad positions to BED format with 10bp windows.
# Creates ±5 bp windows around dyad center for focused analysis.

import sys
import os

if __name__ == "__main__":
    input_file = sys.argv[1]
    file_name = input_file.split('.bed')[0]
    output_file = file_name + "_10bp.bed"
    
    print(f"Output: {output_file}")
    
    with open(input_file, 'r') as ih, open(output_file, 'w') as oh:
        for line in ih.readlines():
            lst = list(line.strip().split('\t'))
            chr_val = lst[0]
            dyad = int(lst[1])
            entry = "{0}\t{1}\t{2}\n".format(chr_val, dyad-5, dyad+5)
            oh.write(entry)
