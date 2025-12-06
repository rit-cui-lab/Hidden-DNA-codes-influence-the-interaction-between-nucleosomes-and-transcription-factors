#!/usr/bin/env python3
# generate_147bp_bed.py
# Author: Chris Carson
# Description: Converts dyad positions to BED format with 147bp windows.
# Creates ±73-74 bp windows around dyad center (nucleosome coordinates).

import sys

if __name__ == "__main__":
    input_file = sys.argv[1]
    file_name = input_file.split('.')[0]
    output_file = file_name + "_147bp.bed"
    
    with open(input_file, 'r') as ih, open(output_file, 'w') as oh:
        for line in ih.readlines():
            lst = list(line.strip().split('\t'))
            chr_val = lst[0]
            dyad = int(lst[1])
            entry = "{0}\t{1}\t{2}\n".format(chr_val, dyad-74, dyad+73)
            oh.write(entry)
