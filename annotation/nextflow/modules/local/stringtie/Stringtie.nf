process STRINGTIE {
	label 'process_medium'
	label 'stringtie'

	input:
	val(all_bam_short)
	val(all_bam_long)

	output:
    path("*.gtf")          , emit: gtf
	path("*.abundance.txt"), emit: abundance
	path("stringtie.log")  , emit: log
    path "versions.yml"    , emit: versions
 
	script:
	if ( all_bam_long != "EMPTY" && all_bam_short != "EMPTY") {
		"""
		stringtie ${all_bam_short} ${all_bam_long} --mix -o stringtie.gtf -p ${task.cpus} -v -A stringtie.gene.abundance.txt
		mv .command.log stringtie.log

		cat <<-VERSIONS > versions.yml
		"${task.process}":
			stringtie: \$( stringtie --version )
		VERSIONS
		"""		
	} else if ( all_bam_short != "EMPTY" && all_bam_long == "EMPTY"){
		"""
		stringtie ${all_bam_short} -o stringtie.gtf -p ${task.cpus} -v -A stringtie.gene.abundance.txt
		mv .command.log stringtie.log

		cat <<-VERSIONS > versions.yml
		"${task.process}":
			stringtie: \$( stringtie --version )
		VERSIONS
		"""
	} else if ( all_bam_long != "EMPTY" && all_bam_short == "EMPTY" ){
		"""
		stringtie ${all_bam_long} -L -o stringtie.gtf -p ${task.cpus} -v -A stringtie.gene.abundance.txt
		mv .command.log stringtie.log

		cat <<-VERSIONS > versions.yml
		"${task.process}":
			stringtie: \$( stringtie --version )
		VERSIONS
		"""
	}
}
