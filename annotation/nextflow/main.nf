#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if (params.help) { exit 0, helpMSG() }

include { ANNOTATO } from './workflows/annotato.nf'

/* 
 Terminal prints
*/

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $workflow.start"
println "Workdir location:"
println "  $workflow.workDir"
println "Launchdir location:"
println "  $workflow.launchDir"
println "Configuration files:"
println "  $workflow.configFiles"
println "Cmd line:"
println "  $workflow.commandLine\u001B[0m"
if (workflow.repository != null){ println "\033[2mGit info: $workflow.repository - $workflow.revision [$workflow.commitId]\u001B[0m" }
println " "
if (workflow.profile.contains('standard') || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.max_cpus, maximal CPUs to use: $params.max_cores\u001B[0m"
    println " "
}

def outdir = new File(params.outdir)
if ( outdir.exists() ) { 
    println ""
    println "\033[0;33mWARNING: Output folder exists. Results will be overwritten, you can adjust the output folder using [--outdir]\033[0m\n"
}

def checkTmp = new File(params.tmpdir)
if ( !checkTmp.exists() ) {
  println ""
  println "\033[0;33mWARNING: ${params.tmpdir} does not exists, creating it.\033[0m\n"
  checkTmp.mkdirs()
}

/*
 Run workflow
*/

workflow {
  println "START THE ANALYSIS\n"
	ANNOTATO()
}

def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    ${c_yellow}Usage examples:${c_reset}
    nextflow run main.nf --genome draftfasta.fa[.gz] --species "Abc def" [--rnaseq] [--protein] [...]

    ${c_yellow}Compulsory input:${c_reset}
    ${c_green}--genome${c_reset}                  Draft genome fasta file contain the assembled contigs/scaffolds
    ${c_green}--species${c_reset}                 Species name for the annotation pipeline, e.g. "Drosophila melanogaster"

    ${c_yellow}Optional input:${c_reset}
    ${c_blue}--protein${c_reset}                  Fasta file containing known protein sequences used as an additional information for gene prediction pipeline.
                                                  Ideally this should come from the same species and/or closely related species. [default: $params.protein]
    ${c_blue}--rnaseq${c_reset}                   A CSV file following the pattern: sample_id,R1_path,R2_path,read_type. 
                                                  This could be generated using `gen_input.py`. Run `python gen_input.py --help` for more information. [default: $params.rnaseq]
    ${c_blue}--long_rna_bam${c_reset}             A BAM file for the alignment of long reads (if any) to the draft genome. Noted that the header of the draft genome need
                                                  to be renamed first before alignment otherwise it will causes trouble for AUGUSTUS and funannotate. [default: $params.long_rna_bam]
    ${c_blue}--short_rna_bam${c_reset}            A BAM file for the alignment of short reads (if any) to the draft genome. Noted that the header of the draft genome need
                                                  to be renamed first before alignment otherwise it will causes trouble for AUGUSTUS and funannotate. [default: $params.short_rna_bam]
    ${c_blue}--knownrepeat${c_reset}              Fasta file containing known repeat sequences of the species, this will be used directly for masking (if --skip_denovo_masking)
                                                  or in combination with the denovo masking. [default: $params.knownrepeat]

    ${c_yellow}Output option:${c_reset}
    ${c_blue}--outdir${c_reset}                   Output directory. [default: $params.outdir]
    ${c_blue}--tracedir${c_reset}                 Pipeline information. [default: $params.tracedir]
    ${c_blue}--publish_dir_mode${c_reset}         Option for nextflow to move data to the output directory. [default: $params.publish_dir_mode]
    ${c_blue}--tmpdir${c_reset}                   Database directory. [default: $params.tmpdir]

    ${c_yellow}Funannotate params:${c_reset}
    ${c_dim}--run_funannotate${c_reset}           Whether to use funannotate for gene prediction. [default: $params.run_funannotate]
    ${c_dim}--organism${c_reset}                  Fungal-specific option. Should be change to "fungus" if the annotated organism is fungal. [default: $params.organism]
    ${c_dim}--ploidy${c_reset}                    Set the ploidy for gene prediction, in case of haploid, a cleaning step will be performed by funannotate to remove
                                                  duplicated contigs/scaffold. [default: $params.ploidy]
    ${c_dim}--buscodb${c_reset}                   BUSCO database used for AUGUSTUS training and evaluation. [default: $params.buscodb]
    ${c_dim}--buscoseed${c_reset}                 AUGUSTUS pre-trained species to start BUSCO. Will be override if rnaseq data is provided. [default: $params.buscoseed]

    ${c_yellow}BRAKER params:${c_reset}
    ${c_dim}--run_braker${c_reset}                Whether to use BRAKER for gene prediction. [default: $params.run_braker]

    ${c_yellow}Skipping options:${c_reset}
    ${c_dim}--skip_rename${c_reset}               Skip renaming genome fasta file by funannotate sort. [default: $params.skip_rename]
    ${c_dim}--skip_all_masking${c_reset}          Skip all masking processes, please be sure that your --genome input is soft-masked before triggering this parameter. [default: $params.skip_all_masking]
    ${c_dim}--skip_denovo_masking${c_reset}       Skip denovo masking using RepeatModeler, this option can only be run when --knownrepeat fasta is provided. [default: $params.skip_denovo_masking]
    ${c_dim}--skip_functional_annotation${c_reset}Skip functional annotation step. [default: $params.skip_functional_annotation]
    ${c_dim}--skip_read_preprocessing${c_reset}   Skip RNASeq preprocessing step. [default: $params.skip_read_preprocessing]

    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
    
    ${c_blue}Engines${c_reset} (choose one):
      conda
      mamba
      docker
      singularity
    
    Per default: -profile slurm,singularity is executed. 
    ${c_reset}
    """.stripIndent()
}