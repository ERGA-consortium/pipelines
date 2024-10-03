process FLAGSTAT {
	label 'process_low'
	label 'minimap2'

	input:
	file(bamFile)
	
	output:
	path("genome_flagstat"), emit: flagstat

	script:
	"""
	samtools flagstat -@ ${task.cpus} *.bam -O json > genome_flagstat
	"""
}
