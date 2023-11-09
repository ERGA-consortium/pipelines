process MINIMAP2INDEX {
	label 'process_medium'
	label 'process_high_memory'
	label 'minimap2'

	input:
	path(fasta)

	output:
	path ("genome.mmi") , emit: index
	path "versions.yml" , emit: versions

	script:
	"""
	minimap2 -x splice -d genome.mmi ${fasta}
	
	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    minimap2: \$(minimap2 --version)
	VERSIONS
	"""
}
