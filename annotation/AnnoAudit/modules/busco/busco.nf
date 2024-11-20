process BUSCO {
    label 'process_medium'
    label 'busco'

    input:
    path(protein)
    val(busco_lineage)

    output:
    path ("annotated")                              , emit: annotated
    path ("annotated/short_summary.specific.*.json"), emit: results

    script:
    if (busco_lineage){
        """
        busco -m prot -i ${protein} -o annotated -c ${task.cpus} -l ${busco_lineage}
        """
    } else {
        """
        busco -m prot -i ${protein} -o annotated --auto-lineage -c ${task.cpus}
        """
    }   
}
