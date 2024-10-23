import sys
import pandas as pd

def parse_gff(gff_file):
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gff_file, sep='\t', header=None, names=columns, comment='#')
    return df

def format_gtf_attributes(attributes, feature):
    attr_dict = {}
    for attribute in attributes.split(';'):
        if '=' in attribute:
            key, value = attribute.strip().split('=', 1)
            attr_dict[key] = value

    gtf_attributes = []
    if 'Parent' in attr_dict:
        gtf_attributes.append(f'gene_id "{attr_dict["Parent"]}"')
    if 'ID' in attr_dict:
        if feature == 'mRNA':
            gtf_attributes.append(f'transcript_id "{attr_dict["ID"]}"')
        elif feature in ['exon', 'CDS']:
            gtf_attributes.append(f'{feature.lower()}_id "{attr_dict["ID"]}"')
            gtf_attributes.append(f'transcript_id "{attr_dict["Parent"]}"')
    
    return '; '.join(gtf_attributes) + ';'

def convert_gff_to_gtf(gff_df):
    gtf_lines = []

    for _, row in gff_df.iterrows():
        if row['feature'] == 'mRNA':
            gtf_attributes = format_gtf_attributes(row['attribute'], 'mRNA')
            gtf_line = f'{row.seqname}\t{row.source}\ttranscript\t{row.start}\t{row.end}\t{row.score}\t{row.strand}\t{row.frame}\t{gtf_attributes}'
            gtf_lines.append(gtf_line)
        elif row['feature'] in ['exon', 'CDS']:
            gtf_attributes = format_gtf_attributes(row['attribute'], row['feature'])
            gtf_line = f'{row.seqname}\t{row.source}\t{row.feature.lower()}\t{row.start}\t{row.end}\t{row.score}\t{row.strand}\t{row.frame}\t{gtf_attributes}'
            gtf_lines.append(gtf_line)
    
    return gtf_lines

def write_gtf_file(gtf_lines, output_file):
    with open(output_file, 'w') as f:
        for line in gtf_lines:
            f.write(line + '\n')

def main(gff_file, output_file):
    gff_df = parse_gff(gff_file)
    gtf_lines = convert_gff_to_gtf(gff_df)
    write_gtf_file(gtf_lines, output_file)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python gff_to_gtf.py input.gff output.gtf")
        sys.exit(1)

    gff_file = sys.argv[1]
    output_file = sys.argv[2]
    main(gff_file, output_file)
