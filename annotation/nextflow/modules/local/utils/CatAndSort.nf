process CATANDSORT {
	label 'process_high'
	label 'process_high_memory'
	label 'minimap2'

	input:
	file(bamFile)
	val(outname)       
	
	output:
	path("${outname}"), emit: all_bam
	path "versions.yml"   , emit: versions

	script:
	"""
	samtools cat -@ ${task.cpus} *.bam | samtools sort -@ ${task.cpus} -o ${outname}

	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	VERSIONS
	"""
}
