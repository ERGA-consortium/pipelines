process DOWNLOAD_OMA {
    label 'process_single'    
     
    output:
    path("LUCA.h5"), emit: database

    script:
    """
    wget https://omabrowser.org/All/LUCA.h5
    """
}
