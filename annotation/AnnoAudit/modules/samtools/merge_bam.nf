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
	bam_count=\$(ls *.bam | wc -l)
	if [ \$bam_count -gt 1 ]; then
		samtools merge -@ ${task.cpus} -o all.sorted.bam *.bam
	else
		mv *.bam all.sorted.bam
	fi 
	"""
}
