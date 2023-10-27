process FASTQC {
	label 'process_low'
	label 'fastqc'
    tag "$meta.sample"

	input:
	tuple val(meta), path(reads)

	output:
    tuple val(meta), path("*.html"), emit: html
    path("*.zip")                  , emit: zip
    path "versions.yml"            , emit: versions

	script:
	"""
    fastqc --noextract -t ${task.cpus} ${reads}

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
	VERSIONS
	"""
}
