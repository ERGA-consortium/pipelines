process PLOT_OMARK {
    label 'process_single'
    label 'seaborn'

    input:
    path(ch_omark_out)
    
    output:
    path("Assessment_plot.pdf"), emit: omark_pdf
    path("Assessment_plot.png"), emit: omark_png

    script:
    """
    python3 ${projectDir}/bin/plot_omark.py --input ${ch_omark_out}
    """
}
