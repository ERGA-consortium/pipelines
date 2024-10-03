process DIAMOND_MAKEDB {
	label 'process_medium'
	label 'diamond'

	input:
	path(protein)

	output:
    path("database*") , emit: database

	script:
	"""
    diamond makedb --in ${protein} -d database
	"""
}