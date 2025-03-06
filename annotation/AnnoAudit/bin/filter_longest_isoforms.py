#!/usr/bin/env python3
import sys
import pandas as pd

def parse_attributes(attr_str):
    """Parse GFF attribute string into a dictionary with vectorized operations"""
    return (
        attr_str.str.extractall(r'(?P<key>[^=;]+)=(?P<value>[^;]+)')
        .droplevel(1)
        .pivot(columns='key', values='value')
    )

def main(input_gff, output_gff):
    # Load GFF with explicit numeric conversion
    df = pd.read_csv(
        input_gff,
        sep='\t',
        comment='#',
        header=None,
        names=[
            'seqid', 'source', 'type', 'start_pos', 'end_pos',
            'score', 'strand', 'phase', 'attributes'
        ],
        dtype={
            'seqid': 'category',
            'source': 'category',
            'type': 'category',
            'strand': 'category',
            'phase': 'category',
            'attributes': 'string'
        }
    )
    
    # Force numeric conversion for coordinates
    df['start_pos'] = pd.to_numeric(df['start_pos'])
    df['end_pos'] = pd.to_numeric(df['end_pos'])

    # Parse attributes
    attr_df = parse_attributes(df['attributes'])
    df = pd.concat([df, attr_df], axis=1)

    # Get gene entries
    genes = df[df['type'] == 'gene'].copy()
    
    # Process mRNAs to find longest isoforms
    mrnas = df[df['type'] == 'mRNA'].dropna(subset=['start_pos', 'end_pos'])
    mrnas['length'] = mrnas['end_pos'] - mrnas['start_pos'] + 1
    
    # Find longest transcript per gene
    longest_isoforms = (
        mrnas.sort_values('length', ascending=False)
        .groupby('Parent', observed=True)
        .head(1)
        [['ID', 'Parent']]
    )
    
    # Collect all IDs to keep
    keep_ids = set(longest_isoforms['ID'])
    gene_ids = set(longest_isoforms['Parent'])
    
    # Find all child features
    children_mask = df['Parent'].isin(keep_ids)
    
    # Combine filters
    final_filter = (
        df['ID'].isin(gene_ids) |  # Gene entries
        df['ID'].isin(keep_ids) |  # mRNA entries
        children_mask              # Child features
    )

    # Build final output
    filtered = pd.concat([
        genes[genes['ID'].isin(gene_ids)],
        df[final_filter & (df['type'] != 'gene')]
    ]).sort_index()

    # Format and save
    filtered[df.columns[:9]].to_csv(
        output_gff,
        sep='\t',
        header=False,
        index=False,
        na_rep='.',
        encoding='utf-8'
    )

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])