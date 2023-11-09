process MULTIQC {
	label 'process_low'
	label 'multiqc'

	input:
	path(qc_input)
	path(config)

	output:
    path "*multiqc_report.html"    , emit: qc_report
	path "*_data"                  , emit: qc_data
    path "versions.yml"            , emit: versions

	script:
	def multiqc_config = config ? "--config ${config}" : ''

	"""
	multiqc --force ${multiqc_config} .

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
	VERSIONS
	"""
}
