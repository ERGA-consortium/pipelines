include { EXTRACT_TRANSCRIPTOME } from '../modules/utils/extract_transcriptome.nf'
include { HISAT2INDEX }           from '../modules/hisat2/hisat2_index.nf'
include { HISAT2 }                from '../modules/hisat2/hisat2.nf'
include { MERGEBAM }              from '../modules/samtools/merge_bam.nf'
include { INDEXSTATS }            from '../modules/samtools/index_stats.nf'

workflow REMAP_TRANSCRIPTOME {
    take:
        ch_genome
        ch_rnaseq
        ch_gff

    main:
        // Extract transcriptome from the input
        EXTRACT_TRANSCRIPTOME ( ch_gff, ch_genome )
        ch_transcriptome = EXTRACT_TRANSCRIPTOME.out.transcriptome
    
        //// HISAT2 index transcriptome
        HISAT2INDEX ( ch_transcriptome )
        transcriptome_index = HISAT2INDEX.out.hisat2_index

        //// HISAT2 mapping back to transcriptome
        HISAT2 ( ch_rnaseq, transcriptome_index )
        ch_out_bam = HISAT2.out.bam

        //// Collect then merging all the bam files
        MERGEBAM ( ch_out_bam.collect() )
        ch_merged_bam = MERGEBAM.out.all_bam

        //// Calculating statistics from all the bam file
        INDEXSTATS ( ch_merged_bam )
        ch_index_stats = INDEXSTATS.out.index_stats
        
    emit:
        ch_merged_bam
        ch_index_stats
}