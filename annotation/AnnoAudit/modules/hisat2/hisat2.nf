process HISAT2 {
    label 'hisat2'
    label 'process_high'
    tag "$meta.sample"

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    path("*transcriptome.bam"), emit: bam

    script:
    if ( meta.paired_end ) {
        """
        hisat2 -x reference.fa -1 ${reads[0]} -2 ${reads[1]} \\
        -p ${task.cpus} --new-summary --summary-file ${meta.sample}.summary | \\
        samtools view -bS | samtools sort -o ${meta.sample}_transcriptome.bam -T tmp --threads ${task.cpus}
        """
    } else {
        """
        hisat2 -x reference.fa -U ${reads[0]} \\
        -p ${task.cpus} --new-summary --summary-file ${meta.sample}.summary | \\
        samtools view -bS | samtools sort -o ${meta.sample}_transcriptome.bam -T tmp --threads ${task.cpus}
        """
    }
    
}