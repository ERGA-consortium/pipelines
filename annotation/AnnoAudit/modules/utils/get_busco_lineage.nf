process GET_BUSCO_LINEAGE {
    label 'process_single'
    label 'biopython'

    input:
    val(query_email)
    val(taxon_id)
    
    output:
    stdout

    script:
    """
    busco_lineage_output=\$(python3 ${projectDir}/bin/get_busco_lineage.py --email ${query_email} --taxon_id ${taxon_id} --busco_lineage_database ${projectDir}/assets/database/busco_lineages.json)
    echo \${busco_lineage_output}
    """
}
