process BLASTX {
	label 'process_medium'
	label 'blast'

	input:
	path(input_fasta)
    path(blastdb)

	output:
    path("blastx_out")  , emit: blastx_out
    path "versions.yml" , emit: versions

	script:
	"""
    blastx -num_threads ${task.cpus} -db combined_prot.noTEs.fa -query ${input_fasta} -out blastx_out

	cat <<-VERSIONS > versions.yml
	"${task.process}":
		blastx: 2.7.1+ )
	VERSIONS
	"""
}
