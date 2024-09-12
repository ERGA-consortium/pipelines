process CALCULATE_STATISTICS {
    label 'process_single'
    label 'pandas'

    input:
    path(ch_fasta)
    path(ch_gff)
    
    output:
    path("statistics.tsv"), emit: statistics

    script:
    """
    python3 ${projectDir}/bin/annot_report.py ${ch_gff} ${ch_fasta}
    """
}
