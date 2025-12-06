#!/usr/bin/env python3
# pfam.py
# Author: Chris Carson
# Description: Comprehensive UniProt protein family analyzer.
# Maps genes to UniProt IDs, fetches Pfam/InterPro families, analyzes by profile labels.

import requests
import pandas as pd
from collections import Counter
import time
import sys
from io import StringIO

# ============================================
# CONFIGURATION
# ============================================

GENE_BATCH_SIZE = 5
FAMILY_BATCH_SIZE = 100
ORGANISM_TAX_ID = '9606'  # 9606 = Human, 10090 = Mouse, 10116 = Rat
BATCH_DELAY_SECONDS = 2

# Input files - UPDATE THESE
CONSENSUS_FILE = 'mammalconsensus.csv'
OUTPUT_LOG = 'analysis_log.txt'
UNIPROT_RESPONSE = 'raw_uniprot_response.tsv'
ANNOTATED_OUTPUT = 'annotated_transcription_factors.csv'
SUMMARY_OUTPUT = 'SUMMARY_by_label.csv'

# ============================================
# LOGGING SETUP
# ============================================

log_file = open(OUTPUT_LOG, 'w', encoding='utf-8')

def log_print(message):
    """Print message to both log file and console."""
    clean_message = message.replace('[OK]', '✓').replace('[ERROR]', '✗')
    print(clean_message, file=log_file)
    log_file.flush()
    print(clean_message)

log_print("="*70)
log_print("UNIPROT PROTEIN FAMILY ANALYZER")
log_print("="*70)

# ============================================
# HELPER FUNCTIONS
# ============================================

def clean_gene_name(gene):
    """Clean gene symbol for UniProt API."""
    if pd.isna(gene):
        return None
    
    cleaned = str(gene).strip().upper()
    cleaned = cleaned.replace('"', '').replace("'", '')
    
    return cleaned if cleaned else None

def get_pfam_names(pfam_ids):
    """Look up human-readable names for Pfam IDs."""
    log_print(f"\nLooking up names for {len(pfam_ids)} Pfam families...")
    
    pfam_names = {}
    
    for pfam_id in pfam_ids:
        try:
            url = f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}'
            response = requests.get(url, timeout=10)
            
            if response.ok:
                data = response.json()
                name = data['metadata']['name']['name']
                pfam_names[pfam_id] = name
            else:
                pfam_names[pfam_id] = pfam_id
            
            time.sleep(0.1)
            
        except Exception:
            pfam_names[pfam_id] = pfam_id
    
    log_print(f"  [OK] Retrieved {len(pfam_names)} Pfam names")
    return pfam_names

def get_interpro_names(interpro_ids):
    """Look up human-readable names for InterPro IDs."""
    log_print(f"\nLooking up names for {len(interpro_ids)} InterPro families...")
    
    interpro_names = {}
    
    for interpro_id in interpro_ids:
        try:
            url = f'https://www.ebi.ac.uk/interpro/api/entry/interpro/{interpro_id}'
            response = requests.get(url, timeout=10)
            
            if response.ok:
                data = response.json()
                name = data['metadata']['name']['name']
                interpro_names[interpro_id] = name
            else:
                interpro_names[interpro_id] = interpro_id
            
            time.sleep(0.1)
            
        except Exception:
            interpro_names[interpro_id] = interpro_id
    
    log_print(f"  [OK] Retrieved {len(interpro_names)} InterPro names")
    return interpro_names

def select_best_uniprot_entry(entries):
    """Choose best UniProt entry (prefer Swiss-Prot/reviewed)."""
    if not entries:
        return None
    
    reviewed = [e for e in entries 
                if e.get('entryType', '').startswith('UniProtKB reviewed')]
    
    best_entry = reviewed[0] if reviewed else entries[0]
    return best_entry['primaryAccession']

def sanitize_filename(label):
    """Make label names safe for filenames."""
    replacements = {
        '/': '_', '\\': '_', ':': '_', '*': '_',
        '?': '_', '"': '_', '<': '_', '>': '_', '|': '_'
    }
    
    safe = label
    for char, replacement in replacements.items():
        safe = safe.replace(char, replacement)
    
    return safe

# ============================================
# GENE â†' UNIPROT ID MAPPING
# ============================================

def map_genes_to_uniprot(genes, batch_size=GENE_BATCH_SIZE):
    """Map gene symbols to UniProt IDs using UniProt API."""
    log_print(f"\n[STEP 1/3] Mapping {len(genes)} genes to UniProt IDs")
    log_print(f"Batch size: {batch_size} genes")
    
    cleaned_genes = []
    for gene in genes:
        cleaned = clean_gene_name(gene)
        if cleaned:
            cleaned_genes.append(cleaned)
    
    log_print(f"Cleaned: {len(cleaned_genes)} valid gene names")
    
    num_batches = (len(cleaned_genes) + batch_size - 1) // batch_size
    log_print(f"Processing in {num_batches} batches")
    
    all_mappings = {}
    
    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, len(cleaned_genes))
        batch = cleaned_genes[start_idx:end_idx]
        
        if batch_num % 10 == 0 or batch_num == num_batches - 1:
            log_print(f"\nBatch {batch_num + 1}/{num_batches}: {batch}")
        
        try:
            response = requests.post(
                'https://rest.uniprot.org/idmapping/run',
                data={
                    'from': 'Gene_Name',
                    'to': 'UniProtKB',
                    'ids': ','.join(batch),
                    'taxId': ORGANISM_TAX_ID
                },
                timeout=30
            )
            
            if not response.ok:
                continue
            
            job_id = response.json()['jobId']
            status_url = f'https://rest.uniprot.org/idmapping/status/{job_id}'
            
            for attempt in range(30):
                time.sleep(2)
                status_response = requests.get(status_url, timeout=30)
                status = status_response.json()
                
                if 'results' in status:
                    break
            
            if 'results' not in status:
                continue
            
            results_data = status['results']
            
            if isinstance(results_data, str):
                results = requests.get(results_data, timeout=60).json()
            elif isinstance(results_data, list):
                results = {'results': results_data}
            else:
                results = status
            
            gene_to_entries = {}
            for item in results.get('results', []):
                gene = item['from']
                if gene not in gene_to_entries:
                    gene_to_entries[gene] = []
                gene_to_entries[gene].append(item['to'])
            
            for gene, entries in gene_to_entries.items():
                uniprot_id = select_best_uniprot_entry(entries)
                if uniprot_id:
                    all_mappings[gene] = uniprot_id
            
            if batch_num % 10 == 0 or batch_num == num_batches - 1:
                total = len(all_mappings)
                log_print(f"  Mapped: {total}/{len(cleaned_genes)} genes so far")
            
            if batch_num < num_batches - 1:
                time.sleep(BATCH_DELAY_SECONDS)
        
        except Exception as e:
            if batch_num % 10 == 0:
                log_print(f"  ERROR: {e}")
            continue
    
    log_print(f"\n[OK] FINAL: {len(all_mappings)}/{len(cleaned_genes)} genes mapped")
    
    final_mapping = {}
    for original_gene in genes:
        cleaned = clean_gene_name(original_gene)
        if cleaned and cleaned in all_mappings:
            final_mapping[original_gene] = all_mappings[cleaned]
    
    if len(final_mapping) < len(genes):
        failed_count = len(genes) - len(final_mapping)
        log_print(f"[WARNING] {failed_count} genes could not be mapped")
    
    return final_mapping

# ============================================
# UNIPROT ID â†' PROTEIN FAMILIES
# ============================================

def fetch_protein_families(uniprot_ids, batch_size=FAMILY_BATCH_SIZE):
    """Fetch Pfam and InterPro family annotations."""
    log_print(f"\n[STEP 2/3] Fetching protein families")
    log_print(f"{len(uniprot_ids)} proteins, batch size: {batch_size}")
    
    if not uniprot_ids:
        log_print("ERROR: No UniProt IDs provided")
        return pd.DataFrame()
    
    num_batches = (len(uniprot_ids) + batch_size - 1) // batch_size
    log_print(f"Processing in {num_batches} batches")
    
    all_dataframes = []
    
    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, len(uniprot_ids))
        batch_ids = uniprot_ids[start_idx:end_idx]
        
        log_print(f"\nFamily batch {batch_num + 1}/{num_batches} ({len(batch_ids)} proteins)")
        
        try:
            query_parts = [f'(accession:{uid})' for uid in batch_ids]
            query = ' OR '.join(query_parts)
            
            response = requests.get(
                'https://rest.uniprot.org/uniprotkb/stream',
                params={
                    'query': query,
                    'format': 'tsv',
                    'fields': 'accession,gene_names,protein_families,xref_pfam,xref_interpro'
                },
                timeout=180
            )
            
            if not response.ok:
                log_print(f"  ERROR: HTTP {response.status_code}")
                continue
            
            batch_df = pd.read_csv(StringIO(response.text), sep='\t')
            
            log_print(f"  [OK] Retrieved {len(batch_df)} proteins")
            all_dataframes.append(batch_df)
            
            if batch_num < num_batches - 1:
                time.sleep(1)
        
        except Exception as e:
            log_print(f"  ERROR: {e}")
            continue
    
    if not all_dataframes:
        log_print("ERROR: No data retrieved from any batch")
        return pd.DataFrame()
    
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    combined_df.to_csv(UNIPROT_RESPONSE, sep='\t', index=False)
    
    log_print(f"\n[OK] Total: {len(combined_df)} proteins")
    log_print(f"[OK] Saved to: {UNIPROT_RESPONSE}")
    
    return combined_df

# ============================================
# ANALYZE BY LABELS
# ============================================

def analyze_by_labels(csv_path, label_column='Consensus', gene_column='Transcription_Factor'):
    """Main analysis function."""
    log_print(f"\n[STEP 3/3] Analyzing {csv_path}")
    
    df = pd.read_csv(csv_path)
    log_print(f"Loaded: {len(df)} rows")
    
    all_genes = df[gene_column].unique().tolist()
    log_print(f"Unique genes: {len(all_genes)}")
    
    gene_to_uniprot = map_genes_to_uniprot(all_genes)
    
    if not gene_to_uniprot or len(gene_to_uniprot) < len(all_genes) * 0.1:
        log_print("ERROR: Insufficient genes mapped - stopping")
        return
    
    uniprot_ids = list(gene_to_uniprot.values())
    family_df = fetch_protein_families(uniprot_ids)
    
    if family_df.empty:
        log_print("ERROR: No family data retrieved")
        return
    
    log_print("\nMerging data...")
    
    df['UniProt_ID'] = df[gene_column].map(gene_to_uniprot)
    
    uniprot_to_families = {}
    for _, row in family_df.iterrows():
        acc = row.get('Entry', '')
        uniprot_to_families[acc] = {
            'Pfam': row.get('Pfam', ''),
            'InterPro': row.get('InterPro', ''),
            'Protein_families': row.get('Protein families', '')
        }
    
    df['Pfam'] = df['UniProt_ID'].map(
        lambda x: uniprot_to_families.get(x, {}).get('Pfam', '')
    )
    df['InterPro'] = df['UniProt_ID'].map(
        lambda x: uniprot_to_families.get(x, {}).get('InterPro', '')
    )
    df['Protein_Family'] = df['UniProt_ID'].map(
        lambda x: uniprot_to_families.get(x, {}).get('Protein_families', '')
    )
    
    df.to_csv(ANNOTATED_OUTPUT, index=False)
    log_print(f"[OK] Saved: {ANNOTATED_OUTPUT}")
    
    log_print(f"\nAnalyzing by '{label_column}' column...")
    
    df_exploded = df.copy()
    df_exploded[label_column] = df_exploded[label_column].astype(str)
    df_exploded['Labels_List'] = df_exploded[label_column].str.
