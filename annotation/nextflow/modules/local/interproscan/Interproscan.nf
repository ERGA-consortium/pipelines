process INTERPROSCAN {
	label 'process_high'

	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://interpro/interproscan:5.64-96.0' :
        'interpro/interproscan:5.64-96.0' }"

	input:
    path(protein)

    output:
    path("${meta.sample}*.bam")       , emit: bam
    path "versions.yml"               , emit: versions

	script:
	"""
	
	
	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    minimap2: \$(minimap2 --version)
		samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	VERSIONS
	"""
}
