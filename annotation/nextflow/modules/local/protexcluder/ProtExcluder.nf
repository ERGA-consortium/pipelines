process PROTEXCLUDER {
	label 'process_medium'
	label 'protexcluder'

	input:
	path(blastx_results)
    path(modeler_output)

	output:
    path("*.fanoProtFinal") , emit: excluded_db
    path "versions.yml"     , emit: versions

	script:
	"""
    ProtExcluder.pl ${blastx_results} ${modeler_output}

	cat <<-VERSIONS > versions.yml
	"${task.process}":
		protexcluder: 1.2.0 )
	VERSIONS
	"""
}
