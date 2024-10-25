process MAKEBLASTDB {
	label 'process_medium'
	label 'blast'

	input:
	path(protein)

	output:
    path("database*") , emit: database

	script:
	"""
    makeblastdb -in ${protein} -dbtype prot -out database
	"""
}