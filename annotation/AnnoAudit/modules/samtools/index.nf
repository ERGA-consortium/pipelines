process INDEX {
	label 'process_medium'
	label 'minimap2'

	input:
	file(ch_genome_bam)
	
	output:
	path("*.bai"), emit: genome_bam_index

	script:
	"""
	samtools index -@ ${task.cpus} ${ch_genome_bam}
	"""
}
