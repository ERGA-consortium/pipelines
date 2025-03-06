process JUNCTION_EXTRACT {
	label 'process_low'
	label 'regtools'

	input:
	path(genome_bam)
	path(genome_bam_index)

	output:
    path("splice_junctions.bed")  , emit: junction_reads

	script:
	"""
    regtools junctions extract -a 8 -m 50 -o splice_junctions.bed ${genome_bam} -s ${params.stranding}
	"""
}
