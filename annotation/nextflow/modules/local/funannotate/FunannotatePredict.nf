process GENE_PREDICTION_FUNANNOTATE {
	label 'process_high'
	label 'process_high_memory'
	label 'funannotate'

	containerOptions = "--bind /env:/env --bind ${params.tmpdir}:/opt/databases --bind /mnt:/mnt"

	input:
	path(softmaskFasta)
	val(protein)
	val(species)
	val(rnabam)
	val(stringtie)

	output:
	path("output")                              , emit: output
	path("output/predict_results/*.proteins.fa"), emit: protseq
	path("output/predict_results/*.gff3")       , emit: gff
	path "versions.yml"                         , emit: versions

	script:
	def args_prot = protein != "EMPTY" ? "--protein_evidence ${protein}" : ""
	def args_rna = rnabam != "EMPTY" ? "--rna_bam ${rnabam}" : ""
	def args_stringtie = stringtie != "EMPTY" ? "--stringtie ${stringtie}" : ""
    def args_buscodb = params.buscodb ? "--busco_db '${params.buscodb}'" : ""
	def args_buscoseed = params.buscoseed == null ? "" : "--busco_seed_species '${params.buscoseed}'"
	def args_ploidy = params.ploidy ? "--ploidy ${params.ploidy}" : ""

	"""
	augustus_species=\$(echo ${species} | sed 's/ /_/g')
	funannotate setup -b ${params.buscodb}

	funannotate predict -i ${softmaskFasta} \\
	-o output \\
	-s "${species}" \\
	--augustus_species \$augustus_species \\
	--cpus ${task.cpus} \\
	--organism ${params.organism} \\
	${args_buscodb} \\
	${args_buscoseed} \\
	${args_prot} \\
	${args_rna} \\
	${args_stringtie} \\
	${args_ploidy}

	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    funannotate: \$( funannotate -version | sed -e "s/funannotate v//g" )
	VERSIONS
	"""
}
