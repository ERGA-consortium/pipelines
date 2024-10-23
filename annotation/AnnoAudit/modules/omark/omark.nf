process OMARK {
    label 'process_single'
    label 'process_medium_memory'
    label 'omark'

    input:
    path(database)     
    path(omamer)
     
    output:
    path("omark_output"), emit: omark
    path("omark_output/_detailed_summary.txt"), emit: omark_results

    script:
    """
    mkdir omark_output
    omark -f ${omamer} -d ${database} -o omark_output
    """
}
