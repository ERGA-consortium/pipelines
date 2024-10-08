import sys
import pandas as pd

def parse_gff3(gff_file, cds_only=False):
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gff_file, sep='\t', header=None, names=columns)

    def extract_gene_id(row):
        attribute_str = str(row['attribute']).rstrip(';')
        attributes = {}
        for item in attribute_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key] = value
        return attributes.get('ID', '') if row['feature'] == 'mRNA' else attributes.get('Parent', '')

    df['gene_id'] = df.apply(extract_gene_id, axis=1)

    genes = df[df['feature'] == 'mRNA']
    if cds_only:
        exons = df[df['feature'] == 'CDS']
    else:
        exons = df[df['feature'] == 'exon']
    cds = df[df['feature'] == 'CDS']

    return genes, exons, cds

def calculate_statistics(genes, exons, cds, genome_size):
    num_genes = len(genes)

    genes_with_one_exon = exons.groupby('gene_id').size() == 1
    num_genes_with_one_exon = genes_with_one_exon.sum()

    gene_lengths = genes['end'] - genes['start'] + 1
    mean_gene_length = round(gene_lengths.mean(), 2)
    median_gene_length = gene_lengths.median()

    total_exons = len(exons)
    exons_per_gene = exons['gene_id'].value_counts()
    mean_exons_per_gene = round(exons_per_gene.mean(), 2)
    median_exons_per_gene = exons_per_gene.median()

    cds['cds_length'] = cds['end'] - cds['start'] + 1
    total_cds_length = cds.groupby('gene_id')['cds_length'].sum()
    mean_total_cds_length = round(total_cds_length.mean(), 2)
    median_total_cds_length = total_cds_length.median()

    sum_cds_length = cds['cds_length'].sum()
    cds_percentage = round((sum_cds_length / genome_size) * 100, 2)

    exons_per_gene = exons.groupby('gene_id').size()
    num_introns_per_gene = exons_per_gene - 1
    total_introns = num_introns_per_gene.sum()

    exons['exon_length'] = exons['end'] - exons['start'] + 1
    exons['classification'] = exons['exon_length'].apply(lambda x: x % 3)
    exons_length_count = exons['classification'].value_counts()
    exon_3n_proportion = f"{exons_length_count.get(0, 0)} ({exons_length_count.get(0, 0) / total_exons * 100:.2f}%)"
    exon_3n1_proportion = f"{exons_length_count.get(1, 0)} ({exons_length_count.get(1, 0) / total_exons * 100:.2f}%)"
    exon_3n2_proportion = f"{exons_length_count.get(2, 0)} ({exons_length_count.get(2, 0) / total_exons * 100:.2f}%)"

    exons_sorted = exons.sort_values(['gene_id', 'start'])
    exons_sorted['intron_start'] = exons_sorted.groupby('gene_id')['end'].shift() + 1
    exons_sorted['intron_length'] = exons_sorted['start'] - exons_sorted['intron_start']
    exons_sorted = exons_sorted.dropna()

    mean_intron_length = round(exons_sorted['intron_length'].mean(), 2)
    median_intron_length = exons_sorted['intron_length'].median()

    stats = {
        'num_genes': num_genes,
        'num_genes_without_introns': f"{num_genes_with_one_exon} ({round(num_genes_with_one_exon/num_genes*100, 2)}%)",
        'mean_gene_length': mean_gene_length,
        'median_gene_length': median_gene_length,
        'num_exons': total_exons,
        'mean_exons_per_gene': mean_exons_per_gene,
        'median_exons_per_gene': median_exons_per_gene,
        'num_exon_3n': exon_3n_proportion,
        'num_exon_3n1': exon_3n1_proportion,
        'num_exon_3n2': exon_3n2_proportion,
        'mean_cds_length': mean_total_cds_length,
        'median_cds_length': median_total_cds_length,
        'total_cds_length': sum_cds_length,
        'percentage_cds_coverage': f"{cds_percentage}%",
        'num_introns': total_introns,
        'mean_intron_length': mean_intron_length,
        'median_intron_length': median_intron_length
    }

    return stats, exons_sorted

def calculate_genome_size(fasta_file):
    total_length = 0
    genome_seq = ""
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq = line.strip()
                genome_seq += seq
                total_length += len(seq)
    return total_length, genome_seq

def save_sequences(df, genome_seq, output_file, is_intron=False):
    with open(output_file, 'w') as f:
        for _, row in df.iterrows():
            if is_intron:
                start = int(row['intron_start'])
                end = int(row['start']) - 1
                if pd.isna(start) or start < 0:
                    continue
            else:
                start = int(row['start']) - 1
                end = int(row['end'])
            
            seq = genome_seq[start:end] 

            if is_intron:
                f.write(f">{row['gene_id']}_{int(row['intron_start'])}_{row['start']-1}\n")
                f.write(f"{seq}\n")
            else:
                f.write(f">{row['gene_id']}_{row['start']}_{row['end']}\n")
                f.write(f"{seq}\n")
                

def print_stats(stats, output_file):
    with open(output_file, 'w') as f:
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")

def main(gff_file, fasta_file, cds_only):
    genes, exons, cds = parse_gff3(gff_file, cds_only)
    genome_size, genome_seq = calculate_genome_size(fasta_file)
    stats, exons_sorted = calculate_statistics(genes, exons, cds, genome_size)
    print_stats(stats, "statistics.tsv")
    #save_sequences(exons, genome_seq, "exon_sequences.fasta")
    save_sequences(exons_sorted, genome_seq, "intron_sequences.fasta", is_intron=True)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python annot_report_pandas.py gff_file fasta_file CDS_only")
        sys.exit(1)

    gff_file = sys.argv[1]
    fasta_file = sys.argv[2]
    cds_only = sys.argv[3]
    main(gff_file, fasta_file, cds_only)
