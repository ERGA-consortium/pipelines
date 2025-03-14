process GET_BUSCO_LINEAGE {
    label 'process_single'
    label 'biopython'

    input:
    val(query_email)
    val(taxon_id)
    val(odb_version)
    
    output:
    stdout

    script:
    """
    busco_lineage_output=\$(python3 ${projectDir}/bin/get_busco_lineage.py --email ${query_email} --taxon_id ${taxon_id} --busco_lineage_database ${projectDir}/assets/database/lineages_file_versions.tsv --odb_version ${odb_version})
    echo \${busco_lineage_output}
    """
}
