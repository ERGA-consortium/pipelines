process COMBINEFASTA {
	label 'process_low'
    label 'unix'

	input:
	path(fasta1)
	path(fasta2)
    val(ch_type)
	
	output:
	path("combined.fasta*"), emit : combined_fasta

	script:
    if ( ch_type == "uniprot" ) {
        """
        bash ${projectDir}/bin/combine.sh ${fasta1} ${fasta2} uniprot
        """
    } else {
        """
        bash ${projectDir}/bin/combine.sh ${fasta1} ${fasta2} fasta
        """
    }
    
}
