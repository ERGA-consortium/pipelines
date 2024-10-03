process FIND_ORTHOLOG {
    label 'process_single'

    input:
    path(pred_ref_blastp)
    path(ref_pred_blastp)
     
    output:
    path("matching_orthologs.txt"), emit: ortholog

    script:
    """
    bash ${projectDir}/bin/find_matching_orthologs.sh ${pred_ref_blastp} ${ref_pred_blastp}
    """
}
