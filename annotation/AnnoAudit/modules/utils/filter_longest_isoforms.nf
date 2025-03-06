process FILTER_LONGEST_GFF {
    label 'process_single'
    label 'pandas'

    input:
    path(ch_gff)
    
    output:
    path("filtered_longest_isoform.gff"), emit: filtered_gff

    script:
    """
    python3 ${projectDir}/bin/filter_longest_isoforms.py ${ch_gff} filtered_longest_isoform.gff
    """
}
