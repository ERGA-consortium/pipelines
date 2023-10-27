process BUILD_DATABASE {
	label 'process_medium_high'
	label 'repeatmodeler'

	input:
	path(cleanfasta)
	val(species)

	output:
	path("${species}*") , emit: database
	path "versions.yml" , emit: versions

	script:

	"""
	#Build new RepeatModeler BLAST database
	BuildDatabase -name "${species}" ${cleanfasta}
	
	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    BuildDatabase: \$( BuildDatabase -v 1> test; head -1 test | sed -e "s~/opt/RepeatModeler/BuildDatabase - ~~g"; rm test )
	VERSIONS
	"""
}
