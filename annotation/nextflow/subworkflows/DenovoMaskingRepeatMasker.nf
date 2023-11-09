include { BUILD_DATABASE    } from '../modules/local/repeatmasker/BuildDatabase.nf'
include { REPEATMODELER     } from '../modules/local/repeatmasker/RepeatModeler.nf'
include { FILTER_REPEAT     } from '../subworkflows/FilterRepeat.nf'
include { REPEATMASKER      } from '../modules/local/repeatmasker/RepeatMasker.nf'

workflow GENOME_MASKING {
    take:
        ch_genome
        val_species
        ch_protein

    main:
        ch_versions_repeatmasker_processing = Channel.empty()
        fix_species = val_species.replaceAll(' ', '_')

        if ( !params.skip_denovo_masking ) {
            // MODULE: Build Database
            BUILD_DATABASE ( ch_genome, fix_species )
            ch_versions_repeatmasker_processing = ch_versions_repeatmasker_processing.mix(BUILD_DATABASE.out.versions)

            // MODULE: Repeat Modeler
            REPEATMODELER ( BUILD_DATABASE.out.database )
            ch_versions_repeatmasker_processing = ch_versions_repeatmasker_processing.mix(REPEATMODELER.out.versions)

            // MODULE: Filter repeat sequences from known proteins from uniprot
            FILTER_REPEAT ( ch_protein, REPEATMODELER.out.familyFa )
            ch_repeat_lib = FILTER_REPEAT.out.repeat_db
            ch_versions_repeatmasker_processing = ch_versions_repeatmasker_processing.mix(FILTER_REPEAT.out.versions)
        } else {
            ch_repeat_lib = Channel.fromPath("${params.knownrepeat}", checkIfExists: true)
        }
        
        // MODULE: Repeat Masking
        REPEATMASKER ( ch_repeat_lib, ch_genome )
        ch_versions_repeatmasker_processing = ch_versions_repeatmasker_processing.mix(REPEATMASKER.out.versions)

        // OUTPUT
        ch_masked_genome = REPEATMASKER.out.masked_genome

    emit:
        ch_masked_genome
        versions = ch_versions_repeatmasker_processing
}