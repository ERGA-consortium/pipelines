process BUILD_STAR {
	label 'process_medium'
	label 'process_high_memory'
	label 'star'

	input:
	path(fasta)

	output:
	path ("STARindex")  , emit: index
	path "versions.yml" , emit: versions

	script:
	"""
	samtools faidx ${fasta}
	NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai`

	mkdir STARindex
	STAR --runMode genomeGenerate \\
	--genomeFastaFiles ${fasta} \\
	--genomeDir STARindex/ \\
	--runThreadN ${task.cpus} \\
	--genomeSAindexNbases \$NUM_BASES
	
	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    star: \$(STAR --version | sed -e "s/STAR_//g")
		samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
	VERSIONS
	"""
}
