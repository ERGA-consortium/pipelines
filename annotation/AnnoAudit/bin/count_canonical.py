from Bio import SeqIO
import sys

def calculate_splice_stats(input_file, output_file):
    canonical = {'GTAG', 'CTAC', 'GCAG', 'CTGC', 'ATAC', 'GTAT'}
    canonical_count = 0
    non_canonical_count = 0
    
    # Count splice sites
    with open(input_file) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            seq = str(record.seq).upper()
            if seq[:2] + seq[-2:] in canonical:
                canonical_count += 1
            else:
                non_canonical_count += 1
    
    # Calculate total and percentages
    total = canonical_count + non_canonical_count
    canonical_percent = (canonical_count / total * 100) if total > 0 else 0
    non_canonical_percent = (non_canonical_count / total * 100) if total > 0 else 0
    
    # Write results to file
    with open(output_file, 'w') as out:
        out.write(f"num_intron_supported\t{total}\n")
        out.write(f"num_intron_supported_canonical\t{canonical_count} ({canonical_percent:.2f}%)\n")
        out.write(f"num_intron_supported_non_canonical\t{non_canonical_count} ({non_canonical_percent:.2f}%)\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    calculate_splice_stats(input_file, output_file)