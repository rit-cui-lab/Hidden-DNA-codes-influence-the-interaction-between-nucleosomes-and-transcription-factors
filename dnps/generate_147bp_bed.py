#!/usr/bin/env python3
"""
generate_147bp_bed.py
Author: Chris Carson

Summary:
This script converts nucleosome dyad positions to 147bp nucleosome footprints.
Given a BED file with dyad positions (single coordinates), it generates 147bp
regions centered on each dyad (-74bp to +73bp), representing the full nucleosome
footprint. This is part of the dNPS pipeline (part 3) for nucleosome pattern analysis.

Usage:
    python generate_147bp_bed.py <input.bed>

Arguments:
    input.bed - BED file with nucleosome dyad positions (chromosome, dyad_position)

Output:
    Creates <input_basename>_147bp.bed with 147bp regions
"""

import sys

if __name__ == "__main__":
    input = sys.argv[1]
    file_name = input.split('.')[0]
    output = file_name + "_147bp.bed"
    with open(input, 'r') as ih, open(output, 'w') as oh:
        for line in ih.readlines():
            lst = list(line.strip().split('\t'))
            chr = lst[0]
            dyad = int(lst[1])
            entry = "{0}\t{1}\t{2}\n".format(chr, dyad - 74, dyad + 73)
            oh.write(entry)
