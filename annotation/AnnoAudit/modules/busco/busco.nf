process BUSCO {
    label 'process_medium'
    label 'busco'

    input:
    path(protein)

    output:
    path ("annotated")                              , emit: annotated
    path ("annotated/short_summary.specific.*.json"), emit: results

    script:
    if (params.lineage){
        """
        busco -m prot -i ${protein} -o annotated -l ${params.lineage} -c ${task.cpus}
        """
    } else {
        """
        busco -m prot -i ${protein} -o annotated --auto-lineage -c ${task.cpus}
        """
    }   
}
