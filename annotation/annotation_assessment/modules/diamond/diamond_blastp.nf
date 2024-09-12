process DIAMOND_BLASTP {
	label 'process_medium'
	label 'diamond'

	input:
	path(input_fasta)
    path(blastdb)
	val(output_file)

	output:
    path(output_file), emit: diamond_out

	script:
	"""
    diamond blastp --threads ${task.cpus} --db database --query ${input_fasta} --out ${output_file} --max-target-seqs 1 --max-hsps 1 \\
    --outfmt 6 qseqid sseqid qlen slen length pident evalue bitscore --masking 0 --iterate
	"""
}
