process FUNANNOTATE_SORT {
	label 'process_low'
	label 'funannotate'

	input:
	path(draftfasta)
	
	output:
	path("contigs.clean.sort.fa"), emit : fasta
	path("header_lookup_tab.txt"), emit : lookuptab
	path("versions.yml")         , emit : versions

	script:
	"""
	# Sort contigs using awk
	bash ${projectDir}/bin/sortFasta.sh ${draftfasta} contigs.sorted.fa

	# Rerun again but already sorted so only rename
	funannotate sort -i contigs.sorted.fa -o contigs.clean.sort.fa

	# Extract the contigs name before and after renaming
	bash ${projectDir}/bin/extractFastaHeaders.sh contigs.sorted.fa contigs.clean.sort.fa header_lookup_tab.txt

	rm contigs.sorted.fa

	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    funannotate: \$( funannotate -version | sed -e "s/funannotate v//g" )
	VERSIONS	
	"""	

}
