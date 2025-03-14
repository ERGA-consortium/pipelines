process GENERATE_PDF {
    label 'process_single'
    label 'reportlab'

    input:
    path(evaluate_json)
    path(busco_plot)
    path(protein_distribution)
    path(intron_without_stop)
    path(intron_with_stop)
    
    output:
    path("*.pdf")

    script:
    """
    realpath "${busco_plot}" > image_paths.txt
    realpath "${protein_distribution}" >> image_paths.txt
    realpath "${intron_without_stop}" >> image_paths.txt
    realpath "${intron_with_stop}" >> image_paths.txt

    python3 ${projectDir}/bin/generate_PDF.py --json ${evaluate_json} --images image_paths.txt --output ${params.pdf}
    """
}
