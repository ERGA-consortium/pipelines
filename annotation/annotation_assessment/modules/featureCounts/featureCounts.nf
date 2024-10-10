process FEATURECOUNTS {
    label 'process_single'
    label 'subread'

    input:
    path(ch_gtf)
    path(ch_genome_bam)
    
    output:
    path("unsupported_genes.txt"), emit: gene_count

    script:
    """
    featureCounts -T ${task.cpus} -p -a ${ch_gtf} -o gene_counts.txt ${ch_genome_bam}
    awk \'NR>1 {total++; if (\$NF==0) zero++} END {printf \"num_gene_unsupported\\t%d (%.2f%%)\\n\", zero, (zero/total)*100}' gene_counts.txt >> unsupported_genes.txt
    """
}
