process ALIGN_STAR {
    label 'process_high'
    label 'star'
    tag "$meta.sample"

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    path("*.bam")       , emit: bam
    path("*.final.out") , emit: log
    path "versions.yml" , emit: versions

    script:
    if ( meta.paired_end ) {
        """
        STAR --genomeDir ${index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --runThreadN ${task.cpus} \\
        --outSAMstrandField intronMotif \\
        --outSAMtype BAM Unsorted \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.sample} \\
        --outTmpDir "./TMP"

        cat <<-VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        VERSIONS
        """   
    } else {
        """
        STAR --genomeDir ${index} \\
        --readFilesIn ${reads[0]} \\
        --runThreadN ${task.cpus} \\
        --outSAMstrandField intronMotif \\
        --outSAMtype BAM Unsorted \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.sample} \\
        --outTmpDir "./TMP"

        cat <<-VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        VERSIONS
        """
    }
    
}