#!/usr/bin/env python3
# processfactorbook.py
# Author: Chris Carson
# Description: Processes MEME motif files using ENCSR-to-protein mappings.
# Filters motifs by target protein list and optionally by biosample type.

import sys
import os

def extract_encsr_from_motif(motif_line):
    """Extract ENCSR code from MOTIF line"""
    parts = motif_line.split()
    if len(parts) > 1:
        return parts[1].split('_')[0]
    return None

def load_tsv_mapping(tsv_file_path, target_biosample=None):
    """Load TSV mapping file and return dict of ENCSR -> protein_name"""
    encsr_to_protein = {}
    with open(tsv_file_path, 'r') as f:
        header = next(f).strip().split('\t')
        print(f"TSV columns: {header}")
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                biosample = parts[1]
                protein_name = parts[2]
                encsr_code = parts[3]
                
                if target_biosample:
                    if biosample.lower() != target_biosample.lower():
                        continue
                
                encsr_to_protein[encsr_code] = protein_name
    
    return encsr_to_protein

def load_target_proteins(protein_list_path):
    """Load the target proteins list"""
    with open(protein_list_path, 'r') as f:
        return set(line.strip() for line in f)

def get_meme_filename(meme_file_path):
    """Get the base filename without path and extension"""
    base_name = os.path.basename(meme_file_path)
    if base_name.endswith('.meme'):
        base_name = base_name[:-5]
    return base_name

def process_meme_file(meme_file_path, tsv_mapping_path, protein_list_path, output_dir=".", biosample=None):
    """Process MEME file using TSV mapping and protein list"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    meme_name = get_meme_filename(meme_file_path)
    encsr_to_protein = load_tsv_mapping(tsv_mapping_path, biosample)
    
    if biosample:
        print(f"Filtering by biosample: {biosample}")
        print(f"Found {len(encsr_to_protein)} ENCSR entries for {biosample}")
    else:
        print(f"No biosample filter applied")
        print(f"Found {len(encsr_to_protein)} total ENCSR entries")
    
    target_proteins = load_target_proteins(protein_list_path)
    
    found_proteins = set()
    remaining_proteins = target_proteins.copy()
    
    with open(meme_file_path, 'r') as f:
        lines = f.readlines()
    
    header_lines = []
    for i, line in enumerate(lines):
        if line.startswith('MOTIF'):
            header_lines = lines[:i]
            break
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('MOTIF'):
            encsr_code = extract_encsr_from_motif(line)
            
            if encsr_code in encsr_to_protein:
                protein_name = encsr_to_protein[encsr_code]
                
                if protein_name in target_proteins:
                    output_lines = []
                    output_lines.extend(header_lines)
                    output_lines.append(lines[i])
                    
                    i += 1
                    while i < len(lines) and not lines[i].startswith('MOTIF'):
                        output_lines.append(lines[i])
                        i += 1
                    
                    output_filename = os.path.join(output_dir, f"{protein_name}_HUMAN_{meme_name}.meme")
                    
                    with open(output_filename, 'w') as f:
                        f.writelines(output_lines)
                    
                    print(f"Created: {output_filename}")
                    found_proteins.add(protein_name)
                    if protein_name in remaining_proteins:
                        remaining_proteins.remove(protein_name)
                    
                    continue
        
        i += 1
    
    with open(os.path.join(output_dir, 'proteins_found.txt'), 'w') as f:
        for protein in sorted(found_proteins):
            f.write(f"{protein}\n")
    
    with open(os.path.join(output_dir, 'remaining_proteins.txt'), 'w') as f:
        for protein in sorted(remaining_proteins):
            f.write(f"{protein}\n")
    
    print(f"Found {len(found_proteins)} proteins")
    print(f"Missing {len(remaining_proteins)} proteins")

if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 6:
        print("Usage: python script.py <meme_file> <tsv_mapping_file> <protein_list_file> [output_directory] [biosample]")
        sys.exit(1)
    
    meme_file = sys.argv[1]
    tsv_mapping = sys.argv[2]
    protein_list = sys.argv[3]
    output_directory = sys.argv[4] if len(sys.argv) >= 5 else "."
    biosample = sys.argv[5] if len(sys.argv) == 6 else None
    
    print(f"Output directory: {os.path.abspath(output_directory)}")
    if biosample:
        print(f"Filtering by biosample: {biosample}")
    
    process_meme_file(meme_file, tsv_mapping, protein_list, output_directory, biosample)
