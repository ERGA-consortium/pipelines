process MERGEBAM {
	label 'process_medium'
	label 'process_medium_memory'
	label 'minimap2'

	input:
	file(bamFile)
	
	output:
	path("all.sorted.bam"), emit: all_bam

	script:
	"""
	samtools merge -@ ${task.cpus} -o all.sorted.bam *.bam 
	"""
}
