import sys
import pandas as pd

def parse_attributes(attribute_str):
    """Parse GFF attribute string into a dictionary"""
    attr_dict = {}
    for attribute in attribute_str.split(';'):
        attribute = attribute.strip()
        if attribute and '=' in attribute:
            key, value = attribute.split('=', 1)
            attr_dict[key.strip()] = value.strip()
    return attr_dict

def parse_gff(gff_file):
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gff_file, sep='\t', header=None, names=columns, comment='#')
    return df

def convert_gff_to_gtf(gff_df):
    gtf_lines = []
    transcript_to_gene = {}  
    
    for _, row in gff_df.iterrows():
        attr_dict = parse_attributes(row['attribute'])
        feature = row['feature']
        
        if feature == 'gene':
            continue
            
        elif feature == 'mRNA':
            transcript_id = attr_dict.get('ID', '')
            gene_id = attr_dict.get('Parent', '')
            if transcript_id and gene_id:
                transcript_to_gene[transcript_id] = gene_id
            
            attributes = [
                f'gene_id "{gene_id}"',
                f'transcript_id "{transcript_id}"'
            ]
            attr_str = '; '.join(attributes) + ';'
            gtf_line = (
                f"{row.seqname}\t{row.source}\ttranscript\t"
                f"{row.start}\t{row.end}\t{row.score}\t"
                f"{row.strand}\t{row.frame}\t{attr_str}"
            )
            gtf_lines.append(gtf_line)
            
        elif feature in ['exon', 'CDS']:
            parent_transcript = attr_dict.get('Parent', '')
            gene_id = transcript_to_gene.get(parent_transcript, '')
            
            attributes = []
            if gene_id:
                attributes.append(f'gene_id "{gene_id}"')
            
            if feature == 'exon':
                exon_id = attr_dict.get('exon_id', attr_dict.get('ID', ''))
                if exon_id:
                    attributes.append(f'exon_id "{exon_id}"')
            elif feature == 'CDS':
                cds_id = attr_dict.get('ID', '')
                if cds_id:
                    attributes.append(f'cds_id "{cds_id}"')
            
            if parent_transcript:
                attributes.append(f'transcript_id "{parent_transcript}"')
            
            attr_str = '; '.join(attributes) + ';'
            gtf_line = (
                f"{row.seqname}\t{row.source}\t{feature.lower()}\t"
                f"{row.start}\t{row.end}\t{row.score}\t"
                f"{row.strand}\t{row.frame}\t{attr_str}"
            )
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