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

if (params.genome_bam) {
    ch_genome_bam = Channel.fromPath("${params.genome_bam}", checkIfExists: true)
} else if (params.rnaseq) {
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
    error ('Input RNASeq or mapped bam file is not specified!')
}

if (!params.ref_protein && params.query_ncbi_prot) {
    if (!params.ncbi_query_email && !params.taxon_id) {
        error ('Reference protein not provided, querying data from the NCBI Entrez protein database. Required: --ncbi_query_email and --taxon_id!!')
    }
}

/* 
  IMPORT MODULES/SUBWORKFLOWS
*/

// MODULES
include { CALCULATE_STATISTICS } from '../modules/utils/calculate_statistics.nf'
include { EXTRACT_INTRON_STATS } from '../modules/utils/extract_intron_stats.nf'
include { GET_BUSCO_LINEAGE    } from '../modules/utils/get_busco_lineage.nf'
include { BUSCO                } from '../modules/busco/busco.nf'
include { DOWNLOAD_OMA         } from '../modules/utils/download_oma.nf'
include { OMAMER               } from '../modules/omark/omamer.nf'
include { OMARK                } from '../modules/omark/omark.nf'
include { DOWNLOAD_UNI         } from '../modules/utils/download_uni.nf'
include { QUERY_NCBI_PROT      } from '../modules/utils/query_ncbi_prot.nf'
include { EXTRACT_PROTEOME     } from '../modules/utils/extract_proteome.nf'
include { COMBINE_REPORT       } from '../modules/utils/combine_report.nf'
include { COMPARE_DISTRIBUTION } from '../modules/utils/compare_distribution.nf'
include { FLAGSTAT             } from '../modules/samtools/flagstat.nf'
include { CUSTOM_GFF2GTF       } from '../modules/utils/gff2gtf_custom.nf'
include { FEATURECOUNTS        } from '../modules/featureCounts/featureCounts.nf'

// SUBWORKFLOWS
include { BEST_RECIPROCAL_HIT  } from '../subworkflows/best_reciprocal_hit.nf'
include { REMAP_GENOME         } from '../subworkflows/remap_genome.nf'

/*
  RUN WORKFLOW
*/

workflow ANNOAUDIT {
    
    if (params.protein) {
        ch_protein = Channel.fromPath("${params.protein}", checkIfExists: true)
    } else {
        EXTRACT_PROTEOME ( ch_genome, ch_gff )
        ch_protein = EXTRACT_PROTEOME.out.proteome
    }

    // General statistics
    CALCULATE_STATISTICS ( ch_genome, ch_gff, params.cds_only )
    ch_statistics_out = CALCULATE_STATISTICS.out.statistics
    ch_intron_fasta = CALCULATE_STATISTICS.out.intron_fasta

    EXTRACT_INTRON_STATS ( ch_intron_fasta, params.genetic_code ,ch_statistics_out )
    ch_all_statistics = EXTRACT_INTRON_STATS.out.statistics

    // Ortholog analysis
    if (params.lineage) {
        ch_busco_lineage = params.lineage
    } else if (params.ncbi_query_email && params.taxon_id && !params.lineage ) {
        GET_BUSCO_LINEAGE ( params.ncbi_query_email, params.taxon_id )
        ch_busco_lineage = GET_BUSCO_LINEAGE.out
    } else {
        ch_busco_lineage = null
    }

    BUSCO ( ch_protein, ch_busco_lineage )
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

    // Protein analysis

    if (params.ref_protein) {
        ch_reference_proteome = Channel.fromPath("${params.ref_protein}", checkIfExists: true)
    } else if (params.query_ncbi_prot) {
        QUERY_NCBI_PROT (params.ncbi_query_email, params.taxon_id, params.ncbi_query_count, params.ncbi_query_batch)
        ch_reference_proteome = QUERY_NCBI_PROT.out.database
    } else {
        DOWNLOAD_UNI()
        ch_reference_proteome = DOWNLOAD_UNI.out.database
    }

    BEST_RECIPROCAL_HIT ( ch_protein, ch_reference_proteome ) // Predicted, Reference
    ch_brh_out = BEST_RECIPROCAL_HIT.out.brh_pairs

    //// Protein length distribution comparison
    COMPARE_DISTRIBUTION ( ch_brh_out )
    ch_compare_distribution_out = COMPARE_DISTRIBUTION.out.compare_distribution

    // RNASeq analysis
    if (params.rnaseq && !params.genome_bam) {
        REMAP_GENOME (ch_genome, ch_rnaseq)
        ch_genome_bam = REMAP_GENOME.out.ch_merged_bam
    }

    // Statistic on the RNASeq bam file
    FLAGSTAT ( ch_genome_bam ) 
    ch_genome_stat = FLAGSTAT.out.flagstat

    CUSTOM_GFF2GTF ( ch_gff )
    annotation_gtf = CUSTOM_GFF2GTF.out.gtf

    FEATURECOUNTS ( annotation_gtf, ch_genome_bam )
    featureCounts_stats = FEATURECOUNTS.out.gene_count
    
    // Generate report
    // Combined all the information to a final output file
    COMBINE_REPORT ( ch_all_statistics, ch_busco_short, ch_omark_out, ch_brh_out, ch_compare_distribution_out, ch_genome_stat, featureCounts_stats)
    
}
