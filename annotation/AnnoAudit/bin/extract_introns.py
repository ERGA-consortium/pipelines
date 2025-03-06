import sys

def parse_attributes(attr_str):
    attr = {}
    parts = attr_str.strip().split(';')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if ' ' in part:
            key, val = part.split(' ', 1)
        else:
            key = part
            val = ''
        key = key.strip()
        val = val.strip(' "\t')
        attr[key] = val
    return attr

def merge_exons(exons):
    if not exons:
        return []
    sorted_exons = sorted(exons, key=lambda x: x['start'])
    merged = [sorted_exons[0].copy()]
    for exon in sorted_exons[1:]:
        last = merged[-1]
        if exon['start'] <= last['end']:
            merged[-1]['end'] = max(last['end'], exon['end'])
        else:
            merged.append(exon.copy())
    return merged

def calculate_transcript_length(exons):
    return sum(exon['end'] - exon['start'] + 1 for exon in exons)

def merge_overlapping_introns(introns):
    if not introns:
        return []
    sorted_introns = sorted(introns, key=lambda x: (x[0], x[1], x[2]))
    merged = [list(sorted_introns[0])]
    for intron in sorted_introns[1:]:
        last = merged[-1]
        if intron[0] == last[0] and intron[1] <= last[2]:
            last[2] = max(last[2], intron[2])
        else:
            merged.append(list(intron))
    return merged

def main(gtf_file):
    gene_to_transcripts = {}
    transcript_exons = {}

    # First pass: Extract gene_id and transcript_id from transcript lines
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type != 'transcript':
                continue
            attributes = parse_attributes(fields[8])
            gene_id = attributes.get('gene_id')
            transcript_id = attributes.get('transcript_id')
            if not gene_id or not transcript_id:
                continue
            if gene_id not in gene_to_transcripts:
                gene_to_transcripts[gene_id] = []
            gene_to_transcripts[gene_id].append(transcript_id)

    # Second pass: Process exon lines
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type != 'exon':
                continue
            seqname = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get('transcript_id')
            if not transcript_id:
                continue
            if transcript_id not in transcript_exons:
                transcript_exons[transcript_id] = []
            transcript_exons[transcript_id].append({'start': start, 'end': end, 'chr': seqname})

    # Extract introns for all isoforms and merge duplicates
    all_introns = []
    for tx_id, exons in transcript_exons.items():
        sorted_exons = sorted(exons, key=lambda x: x['start'])
        for i in range(1, len(sorted_exons)):
            prev = sorted_exons[i-1]
            curr = sorted_exons[i]
            intron_start_gtf = prev['end'] + 1
            intron_end_gtf = curr['start'] - 1
            if intron_start_gtf > intron_end_gtf:
                continue
            bed_start = intron_start_gtf - 1
            bed_end = intron_end_gtf
            chrom = prev['chr']
            all_introns.append((chrom, bed_start, bed_end))
    merged_isoforms_introns = merge_overlapping_introns(all_introns)

    # Write isoforms_introns.bed
    with open('isoforms_introns.bed', 'w') as f:
        for entry in merged_isoforms_introns:
            f.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\n")

    # Find the longest isoform for each gene
    transcript_lengths = {}
    for tx_id, exons in transcript_exons.items():
        transcript_lengths[tx_id] = calculate_transcript_length(exons)
    longest_isoforms = {}
    for gene_id, tx_ids in gene_to_transcripts.items():
        longest_tx_id = None
        longest_length = 0
        for tx_id in tx_ids:
            if tx_id not in transcript_lengths:
                continue
            tx_length = transcript_lengths[tx_id]
            if tx_length > longest_length:
                longest_length = tx_length
                longest_tx_id = tx_id
        if longest_tx_id:
            longest_isoforms[gene_id] = longest_tx_id

    # Extract introns for the longest isoforms
    longest_isoform_introns = []
    for tx_id in longest_isoforms.values():
        exons = transcript_exons[tx_id]
        sorted_exons = sorted(exons, key=lambda x: x['start'])
        for i in range(1, len(sorted_exons)):
            prev = sorted_exons[i-1]
            curr = sorted_exons[i]
            intron_start_gtf = prev['end'] + 1
            intron_end_gtf = curr['start'] - 1
            if intron_start_gtf > intron_end_gtf:
                continue
            bed_start = intron_start_gtf - 1
            bed_end = intron_end_gtf
            chrom = prev['chr']
            longest_isoform_introns.append((chrom, bed_start, bed_end))
    merged_longest_introns = merge_overlapping_introns(longest_isoform_introns)

    # Write longest_isoform_introns.bed
    with open('longest_isoform_introns.bed', 'w') as f:
        for entry in merged_longest_introns:
            f.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\n")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python extract_introns.py input.gtf")
        sys.exit(1)
    main(sys.argv[1])