process PLOT_DISTRIBUTION {
    label 'process_single'
    label 'seaborn'

    input:
    path(ch_brh_out)
    
    output:
    path("Protein_distribution_2000.pdf"), emit: distribution_pdf
    path("Protein_distribution_2000.png"), emit: distribution_png
    path("Protein_length.pdf")           , emit: length_pdf
    path("Protein_length.png")           , emit: length_png
    path("Protein_distribution.p*")

    script:
    """
    python3 ${projectDir}/bin/plot_protein_distribution.py --input_file ${ch_brh_out}
    python3 ${projectDir}/bin/plot_protein_length.py --input ${ch_brh_out}
    """
}
