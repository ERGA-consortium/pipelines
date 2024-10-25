process GFF2GTF {
    label 'process_single'
    label 'gffread'

    input:
    path(ch_gff)
    
    output:
    path("annotation.gtf"), emit: gtf

    script:
    """
    gffread ${ch_gff} -T -o annotation.gtf
    """
}
