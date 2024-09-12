import sys
import pandas as pd

def parse_gff3(gff_file):
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gff_file, sep='\t', header=None, names=columns)

    # Extract gene_id from attributes
    # def extract_gene_id(row):
    #     attributes = dict(item.split('=') for item in str(row['attribute']).rstrip(';').split(';'))
    #     return attributes.get('ID', '') if row['feature'] == 'mRNA' else attributes.get('Parent', '')
    def extract_gene_id(row):
        attribute_str = str(row['attribute']).rstrip(';')
        attributes = {}
        for item in attribute_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key] = value
        return attributes.get('ID', '') if row['feature'] == 'mRNA' else attributes.get('Parent', '')
    
    df['gene_id'] = df.apply(extract_gene_id, axis=1)

    # Filter for genes, exons, and CDS
    genes = df[df['feature'] == 'mRNA']
    exons = df[df['feature'] == 'exon']
    cds = df[df['feature'] == 'CDS']

    return genes, exons, cds

def calculate_statistics(genes, exons, cds, genome_size):
    """Calculates statistics from the provided DataFrames.

    Args:
        genes (pandas.DataFrame): DataFrame containing gene information.
        exons (pandas.DataFrame): DataFrame containing exon information.
        cds (pandas.DataFrame): DataFrame containing CDS information.

    Returns:
        dict: A dictionary containing the calculated statistics.
    """

    # Number of genes
    num_genes = len(genes)

    # Calculate number of genes with only one exon
    genes_with_one_exon = exons.groupby('gene_id').size() == 1
    num_genes_with_one_exon = genes_with_one_exon.sum()
    #genes_without_introns = genes[~genes['gene_id'].isin(exons['gene_id'])]
    #num_genes_without_introns = len(genes_without_introns)

    # Mean and median gene length
    gene_lengths = genes['end'] - genes['start'] + 1 # Including the end position
    mean_gene_length = round(gene_lengths.mean(), 2)
    median_gene_length = gene_lengths.median()

    # Mean and median number of exons per gene
    total_exons = len(exons)
    exons_per_gene = exons['gene_id'].value_counts()
    mean_exons_per_gene = round(exons_per_gene.mean(), 2)
    median_exons_per_gene = exons_per_gene.median()

    # Mean and median CDS size
    cds['cds_length'] = cds['end'] - cds['start'] + 1 # Including the end position
    total_cds_length = cds.groupby('gene_id')['cds_length'].sum()
    mean_total_cds_length = round(total_cds_length.mean(), 2)
    median_total_cds_length = total_cds_length.median()

    # Calculate total CDS size
    sum_cds_length = cds['cds_length'].sum()
    cds_percentage = round(( sum_cds_length / genome_size ) * 100, 2)

    # Calculate number of introns per gene and the total number of introns
    # Need to check the length of exon 
    exons_per_gene = exons.groupby('gene_id').size()
    num_introns_per_gene = exons_per_gene - 1
    total_introns = num_introns_per_gene.sum()

    # Mean and median intron length (assuming one intron per exon pair)
    # intron_lengths = exons.groupby('gene_id')['start'].diff().shift(-1)
    # intron_lengths = intron_lengths.dropna()
    # mean_intron_length = intron_lengths.mean()
    # median_intron_length = intron_lengths.median()

    # Calculate intron length
    exons_sorted = exons.sort_values(['gene_id', 'start'])
    exons_sorted['intron_start'] = exons_sorted.groupby('gene_id')['end'].shift()
    exons_sorted['intron_length'] = exons_sorted['start'] - exons_sorted['intron_start'] - 1 # Not including the exon start base
    exons_sorted = exons_sorted.dropna()

    # Calculate mean and median intron length
    mean_intron_length = round(exons_sorted['intron_length'].mean(), 2)
    median_intron_length = exons_sorted['intron_length'].median()

    # Calculate percentage of 3n, 3n+1, 3n+2 short intron + long intron
    exons_sorted['classification'] = exons_sorted['intron_length'].apply(lambda x: x % 3)
    exons_sorted['length_category'] = exons_sorted['intron_length'].apply(lambda x: 'short' if x < 120 else 'long')
    short_introns = exons_sorted[exons_sorted['length_category'] == 'short']
    if not short_introns.empty:
        short_counts = short_introns['classification'].value_counts()
        short_total = len(short_introns)
        short_3n_proportion = f"{short_counts.get(0, 0)} ({short_counts.get(0, 0) / short_total * 100:.2f}%)"
        short_3n1_proportion = f"{short_counts.get(1, 0)} ({short_counts.get(1, 0) / short_total * 100:.2f}%)"
        short_3n2_proportion = f"{short_counts.get(2, 0)} ({short_counts.get(2, 0) / short_total * 100:.2f}%)"
    long_introns = exons_sorted[exons_sorted['length_category'] == 'long']
    if not long_introns.empty:
        long_counts = long_introns['classification'].value_counts()
        long_total = len(long_introns)
        long_3n_proportion = f"{long_counts.get(0, 0)} ({long_counts.get(0, 0) / long_total * 100:.2f}%)"
        long_3n1_proportion = f"{long_counts.get(1, 0)} ({long_counts.get(1, 0) / long_total * 100:.2f}%)"
        long_3n2_proportion = f"{long_counts.get(2, 0)} ({long_counts.get(2, 0) / long_total * 100:.2f}%)"
    
    stats = {
        'num_genes': num_genes,
        'num_genes_without_introns': str(num_genes_with_one_exon) + " (" + str(round(num_genes_with_one_exon/num_genes*100, 2) ) + "%)",
        'mean_gene_length': mean_gene_length,
        'median_gene_length': median_gene_length,
        'num_exons': total_exons,
        'mean_exons_per_gene': mean_exons_per_gene,
        'median_exons_per_gene': median_exons_per_gene,
        'mean_cds_length': mean_total_cds_length,
        'median_cds_length': median_total_cds_length,
        'total_cds_length': sum_cds_length,
        'percentage_cds_coverage': str(cds_percentage) + "%",
        'num_introns': total_introns,
        'mean_intron_length': mean_intron_length,
        'median_intron_length': median_intron_length,
        'num_short_intron_<120_3n': short_3n_proportion,
        'num_short_intron_<120_3n1': short_3n1_proportion,
        'num_short_intron_<120_3n2': short_3n2_proportion,
        'num_long_intron_>120_3n': long_3n_proportion,
        'num_long_intron_>120_3n1': long_3n1_proportion,
        'num_long_intron_>120_3n2': long_3n2_proportion
    }

    return stats

def calculate_genome_size(fasta_file):
    """Calculates the genome size from a FASTA file."""

    total_length = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                total_length += len(line.strip())
    return total_length

def print_stats(stats, output_file):
     # Open the output file for writing
    with open(output_file, 'w') as f:
        # Print values
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")

def main(gff_file, fasta_file):
    genes, exons, cds = parse_gff3(gff_file)
    genome_size = calculate_genome_size(fasta_file)
    stats = calculate_statistics(genes, exons, cds, genome_size)
    print_stats(stats, "statistics.tsv")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python annot_report.py gff_file fasta_file")
        sys.exit(1)

    gff_file = sys.argv[1]
    fasta_file = sys.argv[2]
    main(gff_file, fasta_file)