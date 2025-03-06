process PSAURON {
	label 'process_low'
	label 'psauron'

	input:
	path(input_fasta)

	output:
    path("psauron_output.csv")  , emit: psauron_out

	script:
	"""
    psauron -i ${input_fasta} -o psauron_output.csv -p -c
	"""
}
