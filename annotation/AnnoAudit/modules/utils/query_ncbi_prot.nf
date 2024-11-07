process QUERY_NCBI_PROT {
    label 'process_single'
    label 'biopython'  

    input:
    val(email)
    val(taxon_id)
    val(target_count)
    val(batch_size)
     
    output:
    path("protein_sequences.fasta"), emit: database

    script:
    """
    python3 ${projectDir}/bin/query_protein_parallel.py --email ${email} --taxon_id ${taxon_id} --target_count ${target_count}
    """
}
