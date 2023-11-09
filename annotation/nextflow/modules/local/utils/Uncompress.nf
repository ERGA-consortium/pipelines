process UNCOMPRESS {
	label 'process_low'
	label 'unix'

	input:
	path(fasta)
	val(type)
	
	output:
	path("*.fasta"), emit : fasta

	script:
	if ( type == "genome" ) {
		"""
		gzip -dc $fasta > genome.fasta
		"""
	} else if ( type == "protein" ) {
		"""
		gzip -dc $fasta > protein.fasta
		"""
	}
}
