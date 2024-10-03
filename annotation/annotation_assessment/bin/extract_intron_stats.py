import sys
from Bio.Seq import Seq
from Bio import SeqIO

def translate_dna(dna_seq, genetic_code):
    seq = Seq(dna_seq)
    protein_seq = seq.translate(table=genetic_code)
    return str(protein_seq)

def count_stop_codons(input_file, genetic_code):
    total_introns = 0
    short_intron_counts_total = {0: 0, 1: 0, 2: 0}  
    short_intron_counts_with_stop = {0: 0, 1: 0, 2: 0} 
    long_intron_counts_total = {0: 0, 1: 0, 2: 0}
    long_intron_counts_with_stop = {0: 0, 1: 0, 2: 0} 

    for record in SeqIO.parse(input_file, "fasta"):
        length_mod_3 = (len(record.seq) % 3)
        total_introns += 1
        protein_seq = translate_dna(str(record.seq), genetic_code)
        if len(record.seq) < 120:
            short_intron_counts_total[length_mod_3] += 1
            if '*' in protein_seq:
                short_intron_counts_with_stop[length_mod_3] += 1
        else:
            long_intron_counts_total[length_mod_3] += 1
            if '*' in protein_seq:
                long_intron_counts_with_stop[length_mod_3] += 1

    return {
        'total_introns': total_introns,
        'short_intron_counts_total': short_intron_counts_total,
        'short_intron_counts_with_stop': short_intron_counts_with_stop,
        'long_intron_counts_total': long_intron_counts_total,
        'long_intron_counts_with_stop': long_intron_counts_with_stop
    }

def main(intron_file, genetic_code):
    intron_stats = count_stop_codons(intron_file, genetic_code)

    with open("stop_codon_statistics.tsv", 'w') as f:
        for key in [0, 1, 2]:
            f.write(f"num_short_intron_<120_3n{key}\t{intron_stats['short_intron_counts_total'][key]} ({intron_stats['short_intron_counts_total'][key]/intron_stats['total_introns']*100:.2f})%\n")
            f.write(f"num_long_intron_>120_3n{key}\t{intron_stats['long_intron_counts_total'][key]} ({intron_stats['long_intron_counts_total'][key]/intron_stats['total_introns']*100:.2f})%\n")
        for key in [0, 1, 2]:
            f.write(f"num_short_intron_<120_3n{key}_with_stop\t{intron_stats['short_intron_counts_with_stop'][key]} ({intron_stats['short_intron_counts_with_stop'][key]/intron_stats['total_introns']*100:.2f})%\n")
            f.write(f"num_long_intron_>120_3n{key}_with_stop\t{intron_stats['long_intron_counts_with_stop'][key]} ({intron_stats['long_intron_counts_with_stop'][key]/intron_stats['total_introns']*100:.2f})%\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python stop_codon_biopython.py intron_file genetic_code")
        sys.exit(1)

    intron_file = sys.argv[1]
    genetic_code = sys.argv[2]
    main(intron_file, genetic_code)

