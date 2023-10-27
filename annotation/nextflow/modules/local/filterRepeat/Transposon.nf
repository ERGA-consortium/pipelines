process TRANSPOSON {
	label 'process_medium'
	label 'transposon'

	input:
	path(protein_chunk)

	output:
    path("*TPSI.allHits") , emit: allHits
    path("*TPSI.topHits") , emit: topHits
    path "versions.yml"   , emit: versions

	script:
	"""
	ln -s /opt/conda/envs/transposon/share/blast-2.2.26/data/BLOSUM* ./
    transposonPSI.pl ${protein_chunk} prot

	cat <<-VERSIONS > versions.yml
	"${task.process}":
		transposonPSI: 1.0.0 )
	VERSIONS
	"""
}
