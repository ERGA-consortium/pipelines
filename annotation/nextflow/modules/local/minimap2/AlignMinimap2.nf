process ALIGN_MINIMAP2 {
	label 'process_high'
	label 'minimap2'
    tag "$meta.sample"

	input:
    tuple val(meta), path(reads)
    path(index)

    output:
    path("${meta.sample}*.bam")       , emit: bam
    path "versions.yml"               , emit: versions
    path "${meta.sample}_summary.log" , emit: log

	script:
	"""
	minimap2 -a -t ${task.cpus} ${index} ${reads[0]} 2> ${meta.sample}_summary.log | samtools view -bS | samtools sort -o ${meta.sample}.sorted.bam -T tmp --threads ${task.cpus}
	
	cat <<-VERSIONS > versions.yml
	"${task.process}":
	    minimap2: \$(minimap2 --version)
		samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	VERSIONS
	"""
}
