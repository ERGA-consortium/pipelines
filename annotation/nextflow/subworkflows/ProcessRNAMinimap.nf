include { MINIMAP2INDEX  } from '../modules/local/minimap2/BuildMinimap2.nf'
include { ALIGN_MINIMAP2 } from '../modules/local/minimap2/AlignMinimap2.nf'

workflow PROCESS_RNA_MINIMAP2 {
    take:
        ch_genome
        ch_rnaseq

    main:
        ch_version_process_rna = Channel.empty()

        // Create index for the draft genome
        MINIMAP2INDEX ( ch_genome )
        ch_version_process_rna = ch_version_process_rna.mix(MINIMAP2INDEX.out.versions)

        // Align reads to the draft genome
        ALIGN_MINIMAP2 ( ch_rnaseq, MINIMAP2INDEX.out.index )
        ch_version_process_rna = ch_version_process_rna.mix(ALIGN_MINIMAP2.out.versions)

        ch_out_bam = ALIGN_MINIMAP2.out.bam
        ch_multiqc_files = ALIGN_MINIMAP2.out.log.collect().ifEmpty([])
    
    emit:
        bam           = ch_out_bam
        multiqc_files = ch_multiqc_files
        versions      = ch_version_process_rna

}