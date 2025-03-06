process CALCULATE_CANONICAL {
    label 'process_single'
    label 'biopython'

    input:
    path(ch_intron_fasta)
    
    output:
    path("Canonical_stats.txt"), emit: canonical_stats

    script:
    """
    python3 ${projectDir}/bin/count_canonical.py ${ch_intron_fasta} Canonical_stats.txt
    """
}
