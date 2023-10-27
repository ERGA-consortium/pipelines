process FUNANNOTATE_CLEAN {
	label 'process_low'
	label 'funannotate'

	input:
	path(draftfasta)
	
	output:
	path("contigs.clean.fa"), emit : fasta
	path("versions.yml")    , emit : versions

	script:
	"""
	funannotate clean -i $draftfasta -o contigs.clean.fa

	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    funannotate: \$( funannotate -version | sed -e "s/funannotate v//g" )
	VERSIONS	
	"""	

}
