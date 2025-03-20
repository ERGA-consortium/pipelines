process FEATURECOUNTS {
    label 'process_single'
    label 'subread'

    input:
    path(ch_gtf)
    path(ch_genome_bam)
    
    output:
    path("unsupported_genes.txt"), emit: gene_count

    script:
    if (!params.rnaseq_single) {
        """
        featureCounts -T ${task.cpus} -p -a ${ch_gtf} -o gene_counts.txt ${ch_genome_bam}
        awk \'NR>1 {total++; if (\$NF<=1) zero++} END {printf \"num_gene_unsupported\\t%d (%.2f%%)\\n\", zero, (zero/total)*100}' gene_counts.txt >> unsupported_genes.txt

        featureCounts -T ${task.cpus} -p -g exon_id -B -O -a ${ch_gtf} -o exon_counts.txt ${ch_genome_bam}
        awk \'NR>1 {total++; if (\$NF<=1) zero++} END {printf \"num_exon_unsupported\\t%d (%.2f%%)\\n\", zero, (zero/total)*100}' exon_counts.txt >> unsupported_genes.txt
        """
    } else {
        """
        featureCounts -T ${task.cpus} -a ${ch_gtf} -o gene_counts.txt ${ch_genome_bam}
        awk \'NR>1 {total++; if (\$NF<=1) zero++} END {printf \"num_gene_unsupported\\t%d (%.2f%%)\\n\", zero, (zero/total)*100}' gene_counts.txt >> unsupported_genes.txt

        featureCounts -T ${task.cpus} -g exon_id -a ${ch_gtf} -o exon_counts.txt ${ch_genome_bam}
        awk \'NR>1 {total++; if (\$NF<=1) zero++} END {printf \"num_exon_unsupported\\t%d (%.2f%%)\\n\", zero, (zero/total)*100}' exon_counts.txt >> unsupported_genes.txt
        """
    }
    
}
