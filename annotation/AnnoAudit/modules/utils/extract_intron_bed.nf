process EXTRACT_INTRON_BED {
	label 'process_low'
	label 'biopython'

	input:
	path(annotation_gtf)

	output:
    path("isoforms_introns.bed")  , emit: intron_bed
    path("longest_isoform_introns.bed"), emit: longest_isoform_introns

	script:
	"""
    python3 ${projectDir}/bin/extract_introns.py ${annotation_gtf}
	"""
}
