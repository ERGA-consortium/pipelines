process OMAMER {
    label 'process_single'
    label 'process_medium_memory'
    label 'omark'

    input:
    path(database)
    path(protein)      

    output:
    path("omamer_output"), emit: omamer

    script:
    """
    omamer search --db ${database} --query ${protein} --out omamer_output
    """
}
