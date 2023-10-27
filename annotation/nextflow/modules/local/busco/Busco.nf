process BUSCO {
    label 'process_medium'
    label 'busco'

    input:
    path(protein)      

    output:
    path ("annotated")                    , emit: annotated
    path ("annotated/short_summary.*.txt"), emit: results
    path "versions.yml"                   , emit: versions

    script:
    """
    busco -m prot -i ${protein} -o annotated -l ${params.buscodb} -c ${task.cpus}

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    VERSIONS
    """
}