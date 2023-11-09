process REPEATMODELER {
	label 'process_high'
	label 'repeatmodeler'

	input:
    path(database)

	output:
	path("*-families.fa")  , emit: familyFa
	path("*-families.stk") , emit: familyStk
	path("*-rmod.log")     , emit: log
	path("RM*")            , emit: runDetail
	path "versions.yml"    , emit: versions

	script:
	"""
	#Run RepeatModeler on the constructed database
	RepeatModeler -database ${database} -LTRStruct -threads ${task.cpus}

	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    RepeatModeler: \$( RepeatModeler -version | sed -e "s/RepeatModeler version //g")
	VERSIONS
	"""
}