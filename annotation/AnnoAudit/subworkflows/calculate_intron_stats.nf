include { INDEX               }      from '../modules/samtools/index.nf'
include { JUNCTION_EXTRACT    }      from '../modules/regtools/regtools_extract_junction.nf'
include { EXTRACT_INTRON_BED  }      from '../modules/utils/extract_intron_bed.nf'
include { BEDTOOLS_INTERSECT  }      from '../modules/bedtools/bedtools_intersect.nf'
include { BEDTOOLS_GETFASTA   }      from '../modules/bedtools/bedtools_getfasta.nf'
include { CALCULATE_CANONICAL }      from '../modules/utils/calculate_canonical.nf'

workflow CALCULATE_INTRON_STATS {
    take:
        ch_genome_fasta
        ch_genome_bam
        ch_annotation_gtf

    main:
        INDEX ( ch_genome_bam )
        ch_genome_bam_index = INDEX.out.genome_bam_index

        JUNCTION_EXTRACT ( ch_genome_bam, ch_genome_bam_index )
        ch_splice_junction = JUNCTION_EXTRACT.out.junction_reads

        EXTRACT_INTRON_BED ( ch_annotation_gtf )
        ch_intron_bed = EXTRACT_INTRON_BED.out.intron_bed

        BEDTOOLS_INTERSECT ( ch_intron_bed, ch_splice_junction )
        ch_supported_introns = BEDTOOLS_INTERSECT.out.supported_introns

        BEDTOOLS_GETFASTA ( ch_supported_introns, ch_genome_fasta )
        ch_intron_sequences = BEDTOOLS_GETFASTA.out.intron_sequences

        CALCULATE_CANONICAL ( ch_intron_sequences )
        ch_canonical_stats = CALCULATE_CANONICAL.out.canonical_stats
        
    emit:
        ch_canonical_stats
}