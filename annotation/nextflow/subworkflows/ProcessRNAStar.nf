include { BUILD_STAR           } from '../modules/local/star/BuildSTAR.nf'
include { ALIGN_STAR           } from '../modules/local/star/AlignSTAR.nf'

workflow PROCESS_RNA_STAR {
    take:
        ch_genome
        ch_rnaseq

    main:
        ch_version_process_rna = Channel.empty()

        // Create index for the draft genome
        BUILD_STAR ( ch_genome )
        ch_version_process_rna = ch_version_process_rna.mix(BUILD_STAR.out.versions)

        // Align reads to the draft genome
        ALIGN_STAR ( ch_rnaseq, BUILD_STAR.out.index )
        ch_version_process_rna = ch_version_process_rna.mix(ALIGN_STAR.out.versions)

        ch_out_bam = ALIGN_STAR.out.bam
        ch_multiqc_files = ALIGN_STAR.out.log.collect().ifEmpty([])
    
    emit:
        bam           = ch_out_bam
        multiqc_files = ch_multiqc_files
        versions      = ch_version_process_rna
}