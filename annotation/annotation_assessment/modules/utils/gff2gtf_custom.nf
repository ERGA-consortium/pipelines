process CUSTOM_GFF2GTF {
    label 'process_single'
    label 'pandas'

    input:
    path(ch_gff)
    
    output:
    path("annotation.gtf"), emit: gtf

    script:
    """
    python3 ${projectDir}/bin/gff2gtf.py ${ch_gff} annotation.gtf
    """
}
