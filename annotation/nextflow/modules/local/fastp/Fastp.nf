process FASTP_TRIM {
	label 'process_medium'
    label 'fastp'
    tag "${meta.sample}"

	input:
    tuple val(meta), path(reads)

	output:
    tuple val(meta), path("${meta.sample}*.trimmed.fastq.gz") , emit: trimmed_reads
    path("*.json")                                            , emit: json
    path("*.html")                                            , emit: html
    path("*.log")                                             , emit: log
    path "versions.yml"                                       , emit: versions

	script:
    def prefix = task.ext.prefix ? "" : "${meta.sample}"

    if ( meta.paired_end ) {
        """
        fastp \\
            -i ${reads[0]} \\
            -I ${reads[1]} \\
            -o ${prefix}.R1.trimmed.fastq.gz \\
            -O ${prefix}.R2.trimmed.fastq.gz \\
            --thread ${task.cpus} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            2> ${prefix}.fastp.log

        cat <<-VERSIONS > versions.yml
        "${task.process}":
            fastp: \$( fastp -v | sed -e "s/fastp //g" )
        VERSIONS
        """
    } else {
        """
        fastp \\
            -i ${reads[0]} \\
            -o ${prefix}.trimmed.fastq.gz \\
            --thread ${task.cpus} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            2> ${prefix}.fastp.log

        cat <<-VERSIONS > versions.yml
        "${task.process}":
            fastp: \$( fastp -v | sed -e "s/fastp //g" )
        VERSIONS
        """
    }
}
