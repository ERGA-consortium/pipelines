process INDEXSTATS {
	label 'process_low'
	label 'minimap2'

	input:
	file(bamFile)
	
	output:
	path("transcriptome_idxstats.txt"), emit: index_stats

	script:
	"""
	samtools idxstats -@ ${task.cpus} *.bam > transcriptome_idxstats.txt
	"""
}
