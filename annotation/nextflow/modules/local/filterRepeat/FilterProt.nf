process FILTERPROT {
	label 'process_medium'
	label 'gaas'

	input:
	path(topHits)
    path(ch_protein_ref)

	output:
	path("combined_prot.noTEs.fa"), emit: filtered_prot
	path "versions.yml"           , emit: versions

	script:
	"""
	cat *.TPSI.topHits | awk '{{if(\$0 ~ /^[^\\/\\/.*]/) print \$5}}' | sort -u > combined_prot.TPSI.topHits.accessions.txt
	gaas_fasta_removeSeqFromIDlist.pl \\
	-f ${ch_protein_ref} \\
	-l combined_prot.TPSI.topHits.accessions.txt \\
	-o combined_prot.noTEs.fa

	cat <<-VERSIONS > versions.yml
	"${task.process}":
		gaas: 1.2.0 )
	VERSIONS
	"""
}
