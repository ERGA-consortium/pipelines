process BUSCO {
    label 'process_medium'
    label 'busco'

    input:
    path(protein)
    val(busco_lineage)

    output:
    path ("annotated")                              , emit: annotated
    path ("annotated/short_summary.specific.*.json"), emit: results
    path ("annotated/busco_figure.png")             , emit: plot

    script:
    if (params.busco_database && busco_lineage) {
        """
        busco -m prot -i ${protein} -o annotated -c ${task.cpus} -l ${busco_lineage} --offline --download_path ${params.busco_database}/${params.odb_version}
        generate_plot.py -wd annotated
        """
    } else if (busco_lineage) {
        """
        busco -m prot -i ${protein} -o annotated -c ${task.cpus} -l ${busco_lineage}
        generate_plot.py -wd annotated
        """
    } else {
        """
        busco -m prot -i ${protein} -o annotated --auto-lineage -c ${task.cpus}
        generate_plot.py -wd annotated
        """
    }   
}
