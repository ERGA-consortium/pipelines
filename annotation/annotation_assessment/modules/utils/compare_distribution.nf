process COMPARE_DISTRIBUTION {
    label 'process_single'
    label 'pandas'

    input:
    path(ch_brh_out)
    
    output:
    path("Distribution_output.txt"), emit: compare_distribution

    script:
    """
    python3 ${projectDir}/bin/compare_distribution.py --brh_out ${ch_brh_out}
    """
}
