process DOWNLOAD_OMA {
    label 'process_single'    
     
    output:
    path("LUCA.h5"), emit: database

    script:
    """
    wget -P ${params.tmpdir} https://omabrowser.org/All/LUCA.h5
    ln -s ${params.tmpdir}/LUCA.h5 .
    """
}
