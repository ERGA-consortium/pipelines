process DOWNLOAD_UNI {
    label 'process_single'    
     
    output:
    path("uniprot_sprot.fasta"), emit: database

    script:
    """
    wget -P ${params.tmpdir} https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gzip -d ${params.tmpdir}/uniprot_sprot.fasta.gz
    ln -s ${params.tmpdir}/uniprot_sprot.fasta .
    """
}
