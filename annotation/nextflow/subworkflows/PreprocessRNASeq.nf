include { FASTQC as FASTQCPRE  } from '../modules/local/qc/FastQC.nf'
include { FASTQC as FASTQCPOST } from '../modules/local/qc/FastQC.nf'
include { FASTP_TRIM           } from '../modules/local/fastp/Fastp.nf'

workflow PREPROCESS_RNA {
    take:
        ch_rnaseq

    main:
        ch_version_process_rna = Channel.empty()

        // FastQC data preprocessing
        FASTQCPRE ( ch_rnaseq )
        ch_version_process_rna = ch_version_process_rna.mix(FASTQCPRE.out.versions)

        // Clean reads with fastp
        FASTP_TRIM ( ch_rnaseq )
        ch_version_process_rna = ch_version_process_rna.mix(FASTP_TRIM.out.versions)

        // FastQC data postprocesing
        FASTQCPOST ( FASTP_TRIM.out.trimmed_reads )
        ch_version_process_rna = ch_version_process_rna.mix(FASTQCPOST.out.versions)
    
        ch_trimmed_reads = FASTP_TRIM.out.trimmed_reads

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(FASTQCPRE.out.zip.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTP_TRIM.out.json.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQCPOST.out.zip.collect().ifEmpty([]))

    emit:
        trimmed_reads = ch_trimmed_reads
        multiqc_files = ch_multiqc_files.collect()
        versions      = ch_version_process_rna

}