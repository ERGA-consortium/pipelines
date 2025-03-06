process PLOT_INTRON_PHASE {
    label 'process_single'
    label 'seaborn'

    input:
    path(short_with_stop)
    path(short_without_stop)
    
    output:
    path("Short_intron_with_stop_distribution.pdf")
    path("Short_intron_with_stop_distribution.png"), emit: short_with_stop_png
    path("Short_intron_without_stop_distribution.pdf")
    path("Short_intron_without_stop_distribution.png"), emit: short_without_stop_png

    script:
    """
    python3 ${projectDir}/bin/plot_intron_phase.py --short_with_stop ${short_with_stop} --short_without_stop ${short_without_stop}
    """
}
