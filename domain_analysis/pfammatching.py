#!/usr/bin/env python3
# pfammatching.py
# Author: Chris Carson
# Description: Parses CD-Search results and filters for Pfam domains only.
# Maps query numbers to gene symbols and exports domain data.

import pandas as pd
from collections import defaultdict
import re

def is_pfam(accession):
    """
    Check if accession is from Pfam database.
    Pfam accessions start with 'pfam' or 'PF' followed by numbers.
    """
    accession_upper = accession.upper()
    accession_lower = accession.lower()
    return accession_lower.startswith('pfam') or accession_upper.startswith('PF')

def create_query_mapping(gene_list):
    """
    Map query numbers to gene names from original list.
    Assumes genes were queried in order (Q#1, Q#2, etc).
    """
    return {f"Q#{i+1}": gene for i, gene in enumerate(gene_list)}

def parse_cdsearch_with_mapping(cdsearch_file, query_to_gene):
    """
    Parse CD-Search results using pre-defined query to gene mapping.
    ONLY includes Pfam domains - filters other databases.
    """
    domains_data = []
    total_lines = 0
    pfam_count = 0
    other_count = 0
    
    with open(cdsearch_file, 'r') as f:
        for line in f:
            if line.startswith('Query') or line.strip() == '':
                continue
            
            if line.startswith('Q#'):
                parts = line.strip().split('\t')
                total_lines += 1
                
                query = parts[0].split()[0]
                gene = query_to_gene.get(query, query)
                
                if len(parts) >= 8:
                    accession = parts[7] if len(parts) > 7 else ''
                    
                    # Only keep Pfam domains
                    if is_pfam(accession):
                        pfam_count += 1
                        domains_data.append({
                            'Query': query,
                            'Gene': gene,
                            'Hit_Type': parts[1],
                            'Accession': accession,
                            'Domain': parts[8] if len(parts) > 8 else '',
                            'From': parts[3],
                            'To': parts[4],
                            'E_value': float(parts[5]),
                            'Bitscore': float(parts[6])
                        })
                    else:
                        other_count += 1
    
    print(f"\nFiltering Statistics:")
    print(f"  Total domain hits: {total_lines}")
    print(f"  Pfam domains: {pfam_count}")
    print(f"  Other databases (filtered out): {other_count}")
    
    return pd.DataFrame(domains_data)

# ============================================
# CONFIGURATION
# ============================================

GENES = []

CD_SEARCH_FILE = 'hitdata.txt'
OUTPUT_FILE = 'Simplified_Domainrecord_pf_pfam.csv'

# ============================================
# MAIN EXECUTION
# ============================================

if __name__ == "__main__":
    print("="*60)
    print("PFAM DOMAIN FILTERING")
    print("="*60)
    
    # Create gene mapping
    query_to_gene = create_query_mapping(GENES)
    
    print(f"\nInput configuration:")
    print(f"  Genes: {len(GENES)} total")
    print(f"  CD-Search file: {CD_SEARCH_FILE}")
    
    # Parse CD-Search results
    domains_df = parse_cdsearch_with_mapping(CD_SEARCH_FILE, query_to_gene)
    
    print(f"\nFinal Dataset:")
    print(f"  Total Pfam domains: {len(domains_df)}")
    print(f"  Unique genes with Pfam domains: {domains_df['Gene'].nunique()}")
    print(f"  Unique Pfam families: {domains_df['Accession'].nunique()}")
    
    print("\n" + "="*60)
    print("Sample of domains found:")
    print("="*60)
    print(domains_df.head(20))
    
    print("\n" + "="*60)
    print("Most common Pfam domains:")
    print("="*60)
    print(domains_df['Domain'].value_counts().head(15))
    
    print("\n" + "="*60)
    print("Pfam domains per gene:")
    print("="*60)
    domains_per_gene = domains_df.groupby('Gene')['Domain'].count().sort_values(ascending=False)
    print(domains_per_gene)
    
    # Export results
    domains_df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n[OK] Saved to {OUTPUT_FILE}")
