/*
  DEFINE FUNCTIONS
*/

public static boolean checkIfFileHasExtension(String s, String[] extn) {
        return Arrays.stream(extn).anyMatch(entry -> s.endsWith(entry));
    }

/*
  VALIDATE AND IMPORT INPUT DATA
*/

// Mandatory

if (params.genome) {
    String[] fasta_extensions = [".fasta", ".fna", ".fa"]
    String[] gzip_extensions  = [".fasta.gz", ".fna.gz", ".fa.gz"]
    is_fasta_input = checkIfFileHasExtension( params.genome.toString().toLowerCase(), fasta_extensions )
    is_gzip_input  = checkIfFileHasExtension( params.genome.toString().toLowerCase(), gzip_extensions )

    if (is_fasta_input && !is_gzip_input) {
        ch_genome = Channel.fromPath("${params.genome}", checkIfExists: true)
        ch_genome_gz = "EMPTY"
    } 

    if (is_gzip_input && !is_fasta_input) {
        ch_genome_gz = Channel.fromPath("${params.genome}", checkIfExists: true)
    } else {
        ch_genome_gz = "EMPTY"
    }

    if (!is_gzip_input && !is_fasta_input) {
        error ('Input genome must be in fasta or fasta.gz format!')
    }
} else { 
    error ('Input genome is not specified!')
}

if (params.species) {
    val_species = params.species
} else { 
    error ('Species name is required!')
}

// Optional
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
    
    def f = file("${params.rnaseq}")
    longCount = 0
    shortCount = 0
    f.splitCsv().each {row ->
        columnValue = row[3]
        if (columnValue.endsWith("long")) {
            longCount++
        } else if (columnValue.endsWith("short")) {
            shortCount++
        }
    }
} 

if (params.long_rna_bam) {
    ch_long_rna_bam = Channel.fromPath("${params.long_rna_bam}", checkIfExists: true)
} 

if (params.short_rna_bam) {
    ch_short_rna_bam = Channel.fromPath("${params.short_rna_bam}", checkIfExists: true)
} 

if (params.protein) {
    String[] fasta_extensions = [".fasta", ".fna", ".fa"]
    String[] gzip_extensions  = [".fasta.gz", ".fna.gz", ".fa.gz"]
    is_fasta_protein = checkIfFileHasExtension( params.protein.toString().toLowerCase(), fasta_extensions )
    is_gzip_protein  = checkIfFileHasExtension( params.protein.toString().toLowerCase(), gzip_extensions )

    if (is_fasta_protein && !is_gzip_protein) {
        ch_protein = Channel.fromPath("${params.protein}", checkIfExists: true)
        ch_protein_gz = "EMPTY"
    }

    if (is_gzip_protein && !is_fasta_protein) {
        ch_protein_gz = Channel.fromPath("${params.protein}", checkIfExists: true)
    } else {
        ch_protein_gz = "EMPTY"
    }

    if (!is_fasta_protein && !is_gzip_protein) {
        error ('Input protein database must be in fasta or fasta.gz format!')
    }
} else {
    ch_protein = "EMPTY"
    ch_protein_gz = "EMPTY"
}

if ( params.skip_denovo_masking && !params.knownrepeat ) {
    error ('Can only skip denovo masking when known repeat fasta are provided!')
}

if ( params.profile ) { exit 1, "--profile is WRONG use -profile" }

if ( params.knownrepeat ) {
    String[] fasta_extensions = [".fasta", ".fna", ".fa"]
    is_fasta_protein = checkIfFileHasExtension( params.knownrepeat.toString().toLowerCase(), fasta_extensions )
    if ( !is_fasta_protein ) {
        error ('Please provide known repeat database in fasta format!')
    }
}

/* 
  IMPORT MODULES/SUBWORKFLOWS
*/

// MODULES
include { UNCOMPRESS as UNCOMPRESS_GENOME } from '../modules/local/utils/Uncompress.nf'
include { UNCOMPRESS as UNCOMPRESS_PROT   } from '../modules/local/utils/Uncompress.nf'
include { CATANDSORT as CATANDSORTSHORT   } from '../modules/local/utils/CatAndSort.nf'
include { CATANDSORT as CATANDSORTLONG    } from '../modules/local/utils/CatAndSort.nf'
include { GENE_PREDICTION_FUNANNOTATE     } from '../modules/local/funannotate/FunannotatePredict.nf'
include { GENE_PREDICTION_BRAKER          } from '../modules/local/braker/Braker.nf'
include { FUNANNOTATE_CLEAN               } from '../modules/local/funannotate/FunannotateClean.nf'
include { FUNANNOTATE_SORT                } from '../modules/local/funannotate/FunannotateSort.nf'
include { STRINGTIE                       } from '../modules/local/stringtie/Stringtie.nf'
include { RENAMEGFF                       } from '../modules/local/utils/RenameGff.nf'
include { MERGEBAM                        } from '../modules/local/utils/MergeBam.nf'
include { MULTIQC                         } from '../modules/local/qc/MultiQC.nf'
include { BUSCO                           } from '../modules/local/busco/Busco.nf'

// FROM NF-CORE
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// SUBWORKFLOWS
include { PROCESS_RNA_MINIMAP2            } from '../subworkflows/ProcessRNAMinimap.nf'
include { PROCESS_RNA_STAR                } from '../subworkflows/ProcessRNAStar.nf'
include { PREPROCESS_RNA                  } from '../subworkflows/PreprocessRNASeq.nf'
include { GENOME_MASKING                  } from '../subworkflows/DenovoMaskingRepeatMasker.nf'

/*
  RUN WORKFLOW
*/

workflow ANNOTATO {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Uncompress data if they are in gzip format
    if ( ch_genome_gz != "EMPTY" ) {
        UNCOMPRESS_GENOME ( ch_genome_gz, "genome" )
        ch_genome = UNCOMPRESS_GENOME.out.fasta.collect()
    }

    if ( ch_protein_gz != "EMPTY" ) {
        UNCOMPRESS_PROT ( ch_protein_gz, "protein" )
        ch_protein = UNCOMPRESS_PROT.out.fasta
    }

    // Removing duplicated contigs if the genome is haploid
    if ( params.ploidy == 1 ) {
        FUNANNOTATE_CLEAN ( ch_genome )
        ch_versions_repeatmasker_processing = ch_versions_repeatmasker_processing.mix(FUNANNOTATE_CLEAN.out.versions.first().ifEmpty(null))
        ch_clean_genome = FUNANNOTATE_CLEAN.out.fasta
    } else {
        ch_clean_genome = ch_genome
    }

    // Sort and rename contigs name to be shorter than 16 bases
    if ( !params.skip_rename ) {
        FUNANNOTATE_SORT ( ch_clean_genome )
        ch_versions = ch_versions.mix(FUNANNOTATE_SORT.out.versions.ifEmpty(null))
        ch_genome_sort = FUNANNOTATE_SORT.out.fasta.collect()
    } else {
        ch_genome_sort = ch_clean_genome.collect()
    }

    // Processing RNASeq data
    if ( params.rnaseq ) {
        if ( !params.skip_read_preprocessing) {
            PREPROCESS_RNA ( ch_rnaseq )
            ch_versions = ch_versions.mix(PREPROCESS_RNA.out.versions.ifEmpty(null))
            ch_multiqc_files = PREPROCESS_RNA.out.multiqc_files
            ch_trimmed_reads = PREPROCESS_RNA.out.trimmed_reads
        } else {
            ch_trimmed_reads = ch_rnaseq
        }
        
        ch_trimmed_reads
                .filter{ it[0].type == "short" }
                .set { ch_shortReads }
        ch_trimmed_reads
                .filter{ it[0].type == "long" }
                .set { ch_longReads }

        // Mapping if the input contain short reads
        if ( params.short_rna_bam ) {
            ch_all_bam_short = ch_short_rna_bam
            no_all_bam_short = false
            println "\u001B[32mFound available bam file for short reads, skipping the mapping process of short reads\033[0m"
        }
        else {
            if ( shortCount > 0 ) {
                PROCESS_RNA_STAR ( ch_genome_sort, ch_shortReads )
                ch_versions = ch_versions.mix(PROCESS_RNA_STAR.out.versions.ifEmpty(null))
                ch_bam_short = PROCESS_RNA_STAR.out.bam
                ch_multiqc_files = ch_multiqc_files.mix(PROCESS_RNA_STAR.out.multiqc_files)

                CATANDSORTSHORT ( ch_bam_short.collect(), "all_short.sorted.bam")
                ch_all_bam_short = CATANDSORTSHORT.out.all_bam
                ch_versions = ch_versions.mix(CATANDSORTSHORT.out.versions.ifEmpty(null))
                no_all_bam_short = false
            } else {
                no_all_bam_short = true
                ch_all_bam_short = Channel.empty()
            }
        }

        // Mapping if the input contain long reads
        if ( params.long_rna_bam ) {
            ch_all_bam_long = ch_long_rna_bam
            no_all_bam_long = false
            println "\u001B[32mFound available bam file for long reads, skipping the mapping process of long reads\033[0m"
        } else {
            if ( longCount > 0 ) {
                PROCESS_RNA_MINIMAP2 ( ch_genome_sort, ch_longReads )
                ch_versions = ch_versions.mix(PROCESS_RNA_MINIMAP2.out.versions.ifEmpty(null))
                ch_bam_long = PROCESS_RNA_MINIMAP2.out.bam
                ch_multiqc_files = ch_multiqc_files.mix(PROCESS_RNA_MINIMAP2.out.multiqc_files)

                CATANDSORTLONG ( ch_bam_long.collect(), "all_long.sorted.bam")
                ch_all_bam_long = CATANDSORTLONG.out.all_bam
                ch_versions = ch_versions.mix(CATANDSORTLONG.out.versions.ifEmpty(null))
                no_all_bam_long = false
            } else {
                no_all_bam_long = true
                ch_all_bam_long = Channel.empty()
            }
        }

        stringtie_all_bam_long = no_all_bam_long ? "EMPTY" : ch_all_bam_long
        stringtie_all_bam_short = no_all_bam_short ? "EMPTY" : ch_all_bam_short

        // Stringtie
        STRINGTIE ( stringtie_all_bam_short, stringtie_all_bam_long )
        ch_versions = ch_versions.mix(STRINGTIE.out.versions.ifEmpty(null))
        ch_stringtie_gtf = STRINGTIE.out.gtf
        ch_multiqc_files = ch_multiqc_files.mix(STRINGTIE.out.log)

        // Merge all bam together for downstream analysis
        ch_all_bam = ch_all_bam_short.mix( ch_all_bam_long )
        MERGEBAM ( ch_all_bam.collect() )
        ch_versions = ch_versions.mix(MERGEBAM.out.versions.ifEmpty(null))
        ch_rna_bam = MERGEBAM.out.all_bam  

    } else {
        if ( params.long_rna_bam && params.short_rna_bam ) {
            ch_all_bam = ch_long_rna_bam.mix( ch_short_rna_bam )
            MERGEBAM ( ch_all_bam.collect() )
            ch_versions = ch_versions.mix(MERGEBAM.out.versions.ifEmpty(null))
            ch_rna_bam = MERGEBAM.out.all_bam
        } else if ( params.long_rna_bam && !params.short_rna_bam ) {
            ch_rna_bam = ch_long_rna_bam
        } else if ( !params.long_rna_bam && params.short_rna_bam ) {
            ch_rna_bam = ch_short_rna_bam
        } else {
            ch_rna_bam = "EMPTY"
            ch_stringtie_gtf = "EMPTY"   
        }

        if ( !params.run_braker ) {
            ch_short_rna_bam_stringtie = params.short_rna_bam ? ch_short_rna_bam : "EMPTY"
            ch_long_rna_bam_stringtie  = params.long_rna_bam ? ch_long_rna_bam : "EMPTY"

            STRINGTIE ( ch_short_rna_bam_stringtie, ch_long_rna_bam_stringtie )
            ch_versions = ch_versions.mix(STRINGTIE.out.versions.ifEmpty(null))
            ch_stringtie_gtf = STRINGTIE.out.gtf
        }
    }

    // Soft masking using repeatmasker 
    if ( !params.skip_all_masking ) {
        GENOME_MASKING ( ch_genome_sort, val_species, ch_protein )
        ch_masked_genome = GENOME_MASKING.out.ch_masked_genome
        ch_versions = ch_versions.mix(GENOME_MASKING.out.versions.ifEmpty(null))
    } else {
        ch_masked_genome = ch_genome_sort
    }

    // Running structural prediction on the genome fasta

    if ( params.run_funannotate && !params.run_braker ) {
        GENE_PREDICTION_FUNANNOTATE ( ch_masked_genome, ch_protein, val_species, ch_rna_bam, ch_stringtie_gtf )
        ch_versions = ch_versions.mix(GENE_PREDICTION_FUNANNOTATE.out.versions.ifEmpty(null))
        ch_predict_prot = GENE_PREDICTION_FUNANNOTATE.out.protseq
        ch_output_gff = GENE_PREDICTION_FUNANNOTATE.out.gff
    } else {
        GENE_PREDICTION_BRAKER ( ch_masked_genome, ch_protein, val_species, ch_rna_bam )
        ch_predict_prot = GENE_PREDICTION_BRAKER.out.protseq
        ch_versions = ch_versions.mix(GENE_PREDICTION_BRAKER.out.versions.ifEmpty(null))
        ch_output_gff = GENE_PREDICTION_BRAKER.out.gff
    }

    // Evaluate using BUSCO
    BUSCO ( ch_predict_prot )
    ch_versions = ch_versions.mix(BUSCO.out.versions.ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO.out.results)

    // Rename contig back to the original version
    if ( !params.skip_rename ) {
        RENAMEGFF ( FUNANNOTATE_SORT.out.lookuptab, ch_output_gff )
        ch_versions = ch_versions.mix(RENAMEGFF.out.versions.ifEmpty(null))
    } 

    // Dump software versions
    //CUSTOM_DUMPSOFTWAREVERSIONS ( ch_versions.unique().collectFile(name: 'collated_versions.yml') )
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    // MultiQC
    ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    MULTIQC ( ch_multiqc_files.collect(), ch_multiqc_config )


    // Running functional annotation
    //if ( !params.skip_functional_annotation ) {
    //    
    //}


}

/*
 FINISHED
*/
