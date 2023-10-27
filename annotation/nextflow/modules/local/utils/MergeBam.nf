process MERGEBAM {
	label 'process_high'
	label 'process_high_memory'
	label 'minimap2'

	input:
	file(bamFile)
	
	output:
	path("all.sorted.bam"), emit: all_bam
	path "versions.yml"   , emit: versions

	script:
	"""
	samtools merge -@ ${task.cpus} -o all.sorted.bam *.bam 

	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	VERSIONS
	"""
}
