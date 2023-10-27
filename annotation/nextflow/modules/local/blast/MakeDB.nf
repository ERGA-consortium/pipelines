process MAKEBLASTDB {
	label 'process_medium'
	label 'blast'

	input:
	path(filtered_protein)

	output:
    path("combined_prot.noTEs.fa.p*") , emit: database
    path "versions.yml"               , emit: versions

	script:
	"""
	gzip -dc ${filtered_protein} > combined_prot.noTEs.fa
    makeblastdb -in combined_prot.noTEs.fa -dbtype prot

	cat <<-VERSIONS > versions.yml
	"${task.process}":
		makeblastdb: 2.7.1+ )
	VERSIONS
	"""
}
