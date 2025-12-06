#!/usr/bin/env python3
# processmemefiles.py
# Author: Chris Carson
# Description: Extracts and organizes MEME motifs by protein name.
# Creates individual MEME files for each target protein.

import sys
import os

def extract_protein_name(motif_line):
    """Get text after 'MOTIF' and before first underscore"""
    parts = motif_line.split()
    if len(parts) > 1:
        return parts[1].split('_')[0]
    return None

def get_meme_filename(meme_file_path):
    """Get the base filename without path"""
    base_name = os.path.basename(meme_file_path)
    return base_name.split('_')[0]

def process_meme_file(meme_file_path, protein_list_path):
    """Process MEME file and extract target proteins"""
    
    meme_name = get_meme_filename(meme_file_path)
    
    with open(protein_list_path, 'r') as f:
        proteins = set(line.strip() for line in f)
    
    remaining_proteins = proteins.copy()
    found_proteins = set()
    
    with open(meme_file_path, 'r') as f:
        lines = f.readlines()
    
    header_lines = lines[:8]
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('MOTIF'):
            protein_name = extract_protein_name(line)
            
            if protein_name in proteins:
                output_lines = []
                output_lines.extend(header_lines)
                output_lines.append(lines[i])
                
                i += 1
                while i < len(lines) and not lines[i].startswith('URL'):
                    output_lines.append(lines[i])
                    i += 1
                
                output_filename = f"{protein_name}_HUMAN_{meme_name}.meme"
                
                with open(output_filename, 'w') as f:
                    f.writelines(output_lines)
                
                found_proteins.add(protein_name)
                remaining_proteins.remove(protein_name)
        i += 1
    
    with open('proteins_found.txt', 'w') as f:
        for protein in found_proteins:
            f.write(f"{protein}\n")
    
    with open('remaining_proteins.txt', 'w') as f:
        for protein in remaining_proteins:
            f.write(f"{protein}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <meme_file> <protein_list_file>")
        sys.exit(1)
    
    meme_file = sys.argv[1]
    protein_list = sys.argv[2]
    process_meme_file(meme_file, protein_list)
