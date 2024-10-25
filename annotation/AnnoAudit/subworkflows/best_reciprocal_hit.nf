include { MAKEBLASTDB as MAKEBLASTDB_1 } from '../modules/blast/makedb.nf'
include { MAKEBLASTDB as MAKEBLASTDB_2 } from '../modules/blast/makedb.nf'
include { BLASTP as BLASTP_1           } from '../modules/blast/blastp.nf'
include { BLASTP as BLASTP_2           } from '../modules/blast/blastp.nf'

include { DIAMOND_MAKEDB as MAKEDB_1   } from '../modules/diamond/diamond_makedb.nf'
include { DIAMOND_MAKEDB as MAKEDB_2   } from '../modules/diamond/diamond_makedb.nf'
include { DIAMOND_BLASTP as DIAMOND_1  } from '../modules/diamond/diamond_blastp.nf'
include { DIAMOND_BLASTP as DIAMOND_2  } from '../modules/diamond/diamond_blastp.nf'

include { FIND_ORTHOLOG                } from '../modules/utils/find_matching_orthologs.nf'

workflow BEST_RECIPROCAL_HIT {
    take:
        ch_predicted_proteome
        ch_reference_proteome

    main:
        if (params.run_blast) {
            MAKEBLASTDB_1 ( ch_predicted_proteome )
            ch_predicted_database = MAKEBLASTDB_1.out.database

            MAKEBLASTDB_2 ( ch_reference_proteome )
            ch_reference_database = MAKEBLASTDB_2.out.database

            BLASTP_1 ( ch_predicted_proteome, ch_reference_database, "pred_ref_blastp" )
            pred_ref_blastp = BLASTP_1.out.blastp_out

            BLASTP_2 ( ch_reference_proteome, ch_predicted_database, "ref_pred_blastp" )
            ref_pred_blastp = BLASTP_2.out.blastp_out
        } else {
            MAKEDB_1 ( ch_predicted_proteome )
            ch_predicted_database = MAKEDB_1.out.database

            MAKEDB_2 ( ch_reference_proteome )
            ch_reference_database = MAKEDB_2.out.database

            DIAMOND_1 ( ch_predicted_proteome, ch_reference_database, "pred_ref_blastp" )
            pred_ref_blastp = DIAMOND_1.out.diamond_out

            DIAMOND_2 ( ch_reference_proteome, ch_predicted_database, "ref_pred_blastp" )
            ref_pred_blastp = DIAMOND_2.out.diamond_out
        }

        FIND_ORTHOLOG ( pred_ref_blastp, ref_pred_blastp )
        brh_pairs = FIND_ORTHOLOG.out.ortholog
        
    emit:
        pred_ref_blastp
        ref_pred_blastp
        brh_pairs
}