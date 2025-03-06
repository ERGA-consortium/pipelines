process PLOT_DISTRIBUTION {
    label 'process_single'
    label 'seaborn'

    input:
    path(ch_brh_out)
    
    output:
    path("Protein_distribution.pdf"), emit: distribution_pdf
    path("Protein_distribution.png"), emit: distribution_png

    script:
    """
    python3 ${projectDir}/bin/plot_protein_distribution.py --input_file ${ch_brh_out}
    """
}
