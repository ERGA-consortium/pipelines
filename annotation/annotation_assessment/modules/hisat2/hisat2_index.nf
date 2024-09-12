process HISAT2INDEX {
    label 'hisat2'
    label 'process_medium'

    input:
    path(reference)

    output:
    path("reference.fa*.ht2"), emit: hisat2_index

    script:
    """
    hisat2-build -p ${task.cpus} ${reference} reference.fa
    """
}