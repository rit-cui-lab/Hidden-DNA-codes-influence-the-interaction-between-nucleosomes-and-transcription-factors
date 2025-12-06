#!/usr/bin/env python3
"""
binning_bed.py
Author: Chris Carson

Summary:
This script creates 200 binned BED files from transcription factor binding sites.
For each input genomic region, it generates 200 overlapping 10bp windows spanning
-1000bp to +1000bp around the center of the binding site. This creates the spatial
bins used for downstream nucleosome positioning analysis (dNPS pipeline part 1).

Usage:
    python binning_bed.py <input.bed> <tf_name> <cell_line> <output_dir>

Arguments:
    input.bed   - BED file with TF binding sites
    tf_name     - Name of the transcription factor
    cell_line   - Cell line identifier
    output_dir  - Directory for output binned BED files
"""

import sys
import numpy as np
import math

if __name__ == "__main__":
    input = sys.argv[1]
    tf_name = sys.argv[2]
    cell_line = sys.argv[3]
    output_dir = sys.argv[4]

    file_name = input.split('.')[0]

    for i in range(200):
        output = output_dir + '/' + tf_name + '_' + cell_line + '_' + file_name + '_' + str(i + 1) + '.bed'
        with open(input, 'r') as ih, open(output, 'w') as oh:
            for line in ih.readlines():
                lst = line.strip().split('\t')
                chr = lst[0]
                start = int(lst[1])
                end = int(lst[2])

                center = start + math.ceil((end - start) / 2)
                new_start = center - 1000 + 10 * i
                new_end = new_start + 9

                result = "{0}\t{1}\t{2}\n".format(chr, new_start, new_end)
                oh.write(result)
            ih.close()
            oh.close()
