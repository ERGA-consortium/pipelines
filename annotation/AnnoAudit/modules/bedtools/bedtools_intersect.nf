process BEDTOOLS_INTERSECT {
	label 'process_low'
	label 'bedtools'

	input:
	path(intron_bed)
    path(splice_junction_bed)

	output:
    path("filtered_supported_introns.bed")  , emit: supported_introns

	script:
	"""
    bedtools intersect -a ${intron_bed} -b ${splice_junction_bed} -wa -wb > supported_introns.bed
	awk 'BEGIN {FS="\\t"; OFS="\\t"} {count[\$1 FS \$2 FS \$3]++} END {for (key in count) if (count[key] > 1) print key, count[key]}' supported_introns.bed > filtered_supported_introns.bed
	"""
}
