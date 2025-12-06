#!/usr/bin/env python3
# calculatensm_multithread.py
# Author: Chris Carson
# Description: Calculates nucleosome occupancy using Gaussian kernel smoothing.
# Processes chromosome dyad data in parallel chunks with multiprocessing.

import sys
import numpy as np
from array import array
import multiprocessing as mp

def process_chunk(chunk_data):
    """Process a chunk of positions, including overlap regions"""
    positions, values, chunk_positions, window_size = chunk_data
    results = []
    
    # Convert to numpy arrays for vectorization
    pos_array = np.array(positions)
    val_array = np.array(values)
    
    # Process each position in this chunk
    for pos in chunk_positions:
        # Create mask for relevant positions
        mask = (pos_array >= pos - window_size) & (pos_array <= pos + window_size)
        relevant_pos = pos_array[mask]
        relevant_val = val_array[mask]
        
        # Vectorized score calculation using Gaussian kernel
        if len(relevant_pos) > 0:
            diffs = relevant_pos - pos
            exp_terms = np.exp(-np.power(diffs/20, 2)/2)
            score = np.sum(relevant_val * exp_terms)
            results.append((pos, score))
    
    return results

def calculate_nsm_occupancy(chr_name, positions, values, output_file):
    """Process a single chromosome using parallel chunks"""
    WINDOW_SIZE = 73
    
    # Convert positions to numpy array
    positions = np.array(positions)
    
    # Determine chunk size based on number of positions
    n_cores = mp.cpu_count()
    n_positions = len(positions)
    positions_per_chunk = max(1000, n_positions // n_cores)
    
    # Create chunks with overlap regions
    chunks = []
    for i in range(0, n_positions, positions_per_chunk):
        chunk_positions = positions[i:i + positions_per_chunk]
        
        min_pos = chunk_positions[0]
        max_pos = chunk_positions[-1]
        overlap_mask = (positions >= min_pos - WINDOW_SIZE) & (positions <= max_pos + WINDOW_SIZE)
        
        chunks.append((
            positions[overlap_mask],
            np.array(values)[overlap_mask],
            chunk_positions,
            WINDOW_SIZE
        ))
    
    # Process chunks in parallel
    with mp.Pool(processes=min(n_cores, len(chunks))) as pool:
        all_results = []
        for chunk_results in pool.imap_unordered(process_chunk, chunks):
            all_results.extend(chunk_results)
    
    # Sort results by position and write to file
    all_results.sort(key=lambda x: x[0])
    with open(output_file, 'w') as oh:
        for pos, score in all_results:
            oh.write(f"{chr_name}\t{pos-1}\t{pos}\t\t{score}\n")

def main():
    if len(sys.argv) != 2:
        print("Usage: script.py input_file")
        return
        
    input_file = sys.argv[1]
    file_name = input_file.split('.')[0]
    output_file = f"{file_name}.bg"
    
    # Read chromosome data
    pos_lst = array('i')
    val_lst = array('f')
    chr_name = None
    
    with open(input_file, 'r') as ih:
        for line in ih:
            chr_name, dyad, score = line.strip().split()
            pos_lst.append(int(dyad))
            val_lst.append(float(score))
    
    if pos_lst and val_lst:
        calculate_nsm_occupancy(chr_name, pos_lst, val_lst, output_file)

if __name__ == "__main__":
    main()
