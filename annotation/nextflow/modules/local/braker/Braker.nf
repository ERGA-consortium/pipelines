process GENE_PREDICTION_BRAKER {
	label 'process_high'
	label 'braker'

	containerOptions = "--bind /env:/env --bind ${params.tmpdir}:/opt/databases"

	input:
	path(softmaskFasta)
	val(protein)
	val(species)
	val(rna_bam)

	output:
	path("braker/*")         , emit: output
	path("braker/braker.aa") , emit: protseq
	path "versions.yml"    , emit: versions

	script:
    def args_fungus = params.organism == "fungus" ? "--fungus" : ""
    def args_prot = protein != "EMPTY" ? "--prot_seq=${protein}" : "--prot_seq=uniprot.fasta"
    def args_rna_bam = rna_bam != "EMPTY" ? "--bam=${rna_bam}" : ""

	"""
    gzip -dc ${projectDir}/db/uniprot_sprot.noTEs.fasta.gz > uniprot.fasta
    augustus_species=\$(echo ${species} | sed 's/ /_/g')

    braker.pl --genome=${softmaskFasta}\\
    ${args_rna_bam} \\
    ${args_prot} \\
    --softmasking \\
    --species=\$augustus_species \\
    ${args_fungus} \\
    --useexisting \\
    --gff3 \\
    --threads ${task.cpus}

	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    braker: \$( braker.pl --version | sed -e "s/braker.pl version //g" )
	VERSIONS
	"""
}
