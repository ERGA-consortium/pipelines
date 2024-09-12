/*
  VALIDATE AND IMPORT INPUT DATA (ALL MANDATORY)
*/

public static boolean checkIfFileHasExtension(String s, String[] extn) {
        return Arrays.stream(extn).anyMatch(entry -> s.endsWith(entry));
    }

// Genome
if (params.genome) {
    String[] fasta_extensions = [".fasta", ".fna", ".fa"]
    is_fasta_input = checkIfFileHasExtension( params.genome.toString().toLowerCase(), fasta_extensions )

    if (is_fasta_input) {
        ch_genome = Channel.fromPath("${params.genome}", checkIfExists: true)
    } else {
        error ('Input genome must be a .fasta[.fna,.fa] file (i.e., unziped)')
    }
} else { 
    error ('Input genome is not specified!')
}

// Predicted GFF
if (params.gff) {
   ch_gff = Channel.fromPath("${params.gff}", checkIfExists: true)
} else {
   error ('Input GFF is not specified!')
}

if (params.rnaseq) {
    ch_rnaseq = Channel.fromPath("${params.rnaseq}", checkIfExists: true)
                .splitCsv(header: true, sep: ",")
                .map{sample_info ->
                    def paired_end = sample_info["R2_path"] ? true : false
                    def read1      = file(sample_info["R1_path"], checkIfExists: true)
                    def read2      = paired_end ? file(sample_info["R2_path"] , checkIfExists: true) : ""
                    def meta       = [:]
                        meta.sample     = sample_info["sample_id"]
                        meta.paired_end = paired_end
                        meta.type       = sample_info["read_type"]
                    return meta.paired_end ? [ meta, [read1, read2] ] : [ meta, [read1] ]
                }
} else {
    error ('Input RNASeq is not specified!')
}

/* 
  IMPORT MODULES/SUBWORKFLOWS
*/

// MODULES
include { CALCULATE_STATISTICS } from '../modules/utils/calculate_statistics.nf'
include { BUSCO                } from '../modules/busco/busco.nf'
include { DOWNLOAD_OMA         } from '../modules/utils/download_oma.nf'
include { OMAMER               } from '../modules/omark/omamer.nf'
include { OMARK                } from '../modules/omark/omark.nf'
include { DOWNLOAD_UNI         } from '../modules/utils/download_uni.nf'
include { EXTRACT_PROTEOME     } from '../modules/utils/extract_proteome.nf'
include { COMBINE_REPORT       } from '../modules/utils/combine_report.nf'
include { COMPARE_DISTRIBUTION } from '../modules/utils/compare_distribution.nf'

// SUBWORKFLOWS
include { BEST_RECIPROCAL_HIT  } from '../subworkflows/best_reciprocal_hit.nf'
include { REMAP_TRANSCRIPTOME  } from '../subworkflows/remap_transcriptome.nf'

/*
  RUN WORKFLOW
*/

workflow ASSESS {
    
    if (params.protein) {
        ch_protein = Channel.fromPath("${params.protein}", checkIfExists: true)
    } else {
        EXTRACT_PROTEOME ( ch_genome, ch_gff )
        ch_protein = EXTRACT_PROTEOME.out.proteome
    }

    //////////////////////////
    // Calculate general statistics
    CALCULATE_STATISTICS ( ch_genome, ch_gff )
    ch_statistics_out = CALCULATE_STATISTICS.out.statistics

    // BUSCO
    BUSCO ( ch_protein )
    ch_busco_short = BUSCO.out.results

    // Check database then run OMARK
    if (params.oma_database) {
      ch_oma_database = Channel.fromPath("${params.oma_database}", checkIfExists: true)
    } else {
      DOWNLOAD_OMA()
      ch_oma_database = DOWNLOAD_OMA.out.database
    }

    OMAMER (ch_oma_database, ch_protein)

    OMARK (ch_oma_database, OMAMER.out.omamer)
    ch_omark_out = OMARK.out.omark_results

    //////////////////////////
    // Protein
    //// BRH (with Diamond or Blast)

    if (params.ref_protein) {
        ch_reference_proteome = Channel.fromPath("${params.ref_protein}", checkIfExists: true)
    } else {
        DOWNLOAD_UNI()
        ch_reference_proteome = DOWNLOAD_UNI.out.database
    }

    BEST_RECIPROCAL_HIT ( ch_protein, ch_reference_proteome ) // Predicted, Reference
    ch_brh_out = BEST_RECIPROCAL_HIT.out.brh_pairs

    //// Protein length distribution comparison
    COMPARE_DISTRIBUTION ( ch_brh_out )
    ch_compare_distribution_out = COMPARE_DISTRIBUTION.out.compare_distribution

    //////////////////////////
    // Re-mapping with RNASeq data
    REMAP_TRANSCRIPTOME ( ch_genome, ch_rnaseq, ch_gff )
    ch_transcriptome_out = REMAP_TRANSCRIPTOME.out.ch_index_stats

    //////////////////////////
    // Combined all the information to a final output file
    COMBINE_REPORT ( ch_statistics_out, ch_busco_short, ch_omark_out, ch_brh_out, ch_compare_distribution_out, ch_transcriptome_out )
    
}
