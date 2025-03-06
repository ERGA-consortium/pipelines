process EXTRACT_INTRON_STATS {
    label 'process_single'
    label 'biopython'

    input:
    path(ch_fasta)
    val(genetic_code)
    path(general_stats)
    
    output:
    path("statistics.tsv"), emit: statistics
    path("short_with_stop.pickle"), emit: short_with_stop
    path("short_without_stop.pickle"), emit: short_without_stop

    script:
    """
    python3 ${projectDir}/bin/extract_intron_stats.py ${ch_fasta} ${genetic_code}
    cat stop_codon_statistics.tsv >> statistics.tsv
    """
}
