process DOWNLOAD_UNI {
    label 'process_single'    
     
    output:
    path("uniprot_sprot.fasta"), emit: database

    script:
    """
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gzip -d uniprot_sprot.fasta.gz
    """
}
