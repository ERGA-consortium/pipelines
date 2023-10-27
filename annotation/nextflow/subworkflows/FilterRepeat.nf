include { COMBINEFASTA as COMBINEUNIPROT } from '../modules/local/utils/CombineFasta.nf'
include { COMBINEFASTA as COMBINEKNOWN   } from '../modules/local/utils/CombineFasta.nf'
include { TRANSPOSON                     } from '../modules/local/filterRepeat/Transposon.nf'
include { FILTERPROT                     } from '../modules/local/filterRepeat/FilterProt.nf'
include { MAKEBLASTDB                    } from '../modules/local/blast/MakeDB.nf'
include { BLASTX                         } from '../modules/local/blast/Blastx.nf'
include { PROTEXCLUDER                   } from '../modules/local/protexcluder/ProtExcluder.nf'

workflow FILTER_REPEAT {
    take:
        ch_protein
        ch_familyfa

    main:
        ch_version_process_filter = Channel.empty()
        // Provide already identified transposon data from the Uniprot data base
        ch_uniprot_ref = Channel.fromPath("${projectDir}/db/uniprot_sprot.noTEs.fasta.gz")

        // If there is user provided protein data, we can use it as well
        if ( ch_protein != "EMPTY" ) {
            ch_protein_ref_chunk = ch_protein.splitFasta( by: 200, file: true)

            // First identify sequences that contain transposon
            TRANSPOSON ( ch_protein_ref_chunk )

            // Then filtered those sequences
            FILTERPROT ( TRANSPOSON.out.topHits.collect(), ch_protein )
            ch_version_process_filter = ch_version_process_filter.mix(FILTERPROT.out.versions)

            // Finally, combine the output with the precalculated one from uniprot
            ch_type = "uniprot"
            COMBINEUNIPROT ( ch_uniprot_ref, FILTERPROT.out.filtered_prot, ch_type )

            ch_ref_fasta = COMBINEUNIPROT.out.combined_fasta
        } else {
            ch_ref_fasta = ch_uniprot_ref
        }

        // Make blastdb from the (combined) filtered protein
        MAKEBLASTDB ( ch_ref_fasta )
        ch_version_process_filter = ch_version_process_filter.mix(MAKEBLASTDB.out.versions)

        // Blast the denovo repeat family to the protein database
        BLASTX ( ch_familyfa, MAKEBLASTDB.out.database )
        ch_version_process_filter = ch_version_process_filter.mix(BLASTX.out.versions)

        // Filter out the repeat family that mapped to the protein database
        PROTEXCLUDER ( BLASTX.out.blastx_out, ch_familyfa )
        ch_version_process_filter = ch_version_process_filter.mix(PROTEXCLUDER.out.versions)

        if ( params.knownrepeat ) {
            ch_knownrepeat = Channel.fromPath("${params.knownrepeat}", checkIfExists: true)
            ch_type = "fasta"
            COMBINEKNOWN ( ch_knownrepeat, PROTEXCLUDER.out.excluded_db, ch_type )
            ch_repeat_db = COMBINEKNOWN.out.combined_fasta
        } else {
            ch_repeat_db = PROTEXCLUDER.out.excluded_db
        }

    emit:
        repeat_db = ch_repeat_db
        versions  = ch_version_process_filter

}