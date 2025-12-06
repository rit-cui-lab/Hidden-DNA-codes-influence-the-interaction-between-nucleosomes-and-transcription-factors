#!/usr/bin/env python3
# getproteinseq.py
# Author: Chris Carson
# Description: Downloads protein sequences from NCBI for gene symbols.
# Retrieves sequences in FASTA format and saves to single file.

from Bio import Entrez
import time

# ============================================
# CONFIGURATION
# ============================================

# REQUIRED: Set your email (NCBI policy)
Entrez.email = "your.email@example.com"

GENES = [ ]

OUTPUT_FILE = 'protein_sequences.fasta'
REQUEST_DELAY_SECONDS = 0.4  # Respect NCBI rate limits

# ============================================
# MAIN EXECUTION
# ============================================

if __name__ == "__main__":
    print("="*60)
    print("NCBI PROTEIN SEQUENCE RETRIEVAL")
    print("="*60)
    print(f"\nRetrieving {len(GENES)} protein sequences...")
    print(f"Rate limit: 1 request every {REQUEST_DELAY_SECONDS} seconds")
    
    all_sequences = []
    failed_genes = []
    
    for i, gene in enumerate(GENES):
        # Build search query: gene symbol + human + RefSeq
        search_term = f"{gene}[Gene] AND Homo sapiens[Organism] AND RefSeq[Filter]"
        
        try:
            # Search for protein IDs
            search_handle = Entrez.esearch(db="protein", term=search_term, retmax=1)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            if search_results['IdList']:
                protein_id = search_results['IdList'][0]
                
                # Fetch protein sequence
                fetch_handle = Entrez.efetch(db="protein", id=protein_id, 
                                            rettype="fasta", retmode="text")
                sequence = fetch_handle.read()
                fetch_handle.close()
                
                all_sequences.append(sequence)
                print(f"  [{i+1:3d}/{len(GENES)}] {gene:10s} - Retrieved")
            else:
                failed_genes.append(gene)
                print(f"  [{i+1:3d}/{len(GENES)}] {gene:10s} - No sequence found")
            
            # Respect NCBI rate limits
            time.sleep(REQUEST_DELAY_SECONDS)
            
        except Exception as e:
            print(f"  [{i+1:3d}/{len(GENES)}] {gene:10s} - ERROR: {e}")
            failed_genes.append(gene)
            time.sleep(1)
    
    # ============================================
    # SAVE RESULTS
    # ============================================
    
    with open(OUTPUT_FILE, 'w') as output_file:
        output_file.write(''.join(all_sequences))
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"  Successfully retrieved: {len(all_sequences)}/{len(GENES)}")
    print(f"  Failed: {len(failed_genes)}")
    print(f"  Output file: {OUTPUT_FILE}")
    
    if failed_genes:
        print(f"\nFailed genes:")
        for gene in failed_genes[:10]:
            print(f"  - {gene}")
        if len(failed_genes) > 10:
            print(f"  ... and {len(failed_genes) - 10} more")
    
