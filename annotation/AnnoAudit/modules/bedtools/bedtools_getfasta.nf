process BEDTOOLS_GETFASTA {
	label 'process_low'
	label 'bedtools'

	input:
	path(supported_introns)
    path(genome_fasta)

	output:
    path("intron_sequences.fa")  , emit: intron_sequences

	script:
	"""
    bedtools getfasta -fi ${genome_fasta} -bed ${supported_introns} -fo intron_sequences.fa
	"""
}
