process BLASTP {
	label 'process_high'
	label 'blast'

	input:
	path(input_fasta)
    path(blastdb)
	val(output_file)

	output:
    path(output_file)  , emit: blastp_out

	script:
	"""
    blastp -num_threads ${task.cpus} -db database -query ${input_fasta} -out ${output_file} -max_target_seqs 1 -max_hsps 1 \\
	-outfmt "6 qseqid sseqid qlen slen length pident evalue bitscore"
	"""
}
