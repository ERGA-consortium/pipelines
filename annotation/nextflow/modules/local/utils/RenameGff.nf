process RENAMEGFF {
	label 'process_low'
    label 'process_single'
	label 'tools'

	input:
	file(lookup)
	file(inGff)   
	
	output:
	path("renamed_output.gff"), emit: outGff
	path "versions.yml"       , emit: versions

	script:
	"""
	python ${projectDir}/bin/rename_gff.py -l ${lookup} -i ${inGff} -o renamed_output.gff

	cat <<-VERSIONS > versions.yml
	"${task.process}":
		python: 3.11.6
	    pandas: 2.1.3
	VERSIONS
	"""
}
