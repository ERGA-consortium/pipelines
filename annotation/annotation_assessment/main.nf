#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if (params.help) { exit 0, helpMSG() }

include { ASSESS } from './workflows/assessment.nf'

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
    ASSESS()
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

    ${c_yellow}Input parameter:${c_reset}
    ${c_green}--genome${c_reset}                  Draft genome fasta file contain the assembled contigs/scaffolds
    ${c_green}--gff${c_reset}                     Annotation file that need to be evaluated
    ${c_green}--genome_bam${c_reset}              BAM file contain the mapped information from the RNASeq to the genome FASTA.
    ${c_green}--rnaseq${c_reset}                  A CSV file following the pattern: sample_id,R1_path,R2_path,read_type. Required if `genome_bam` is not provided.


    ${c_yellow}Optional input:${c_reset}
    ${c_blue}--protein${c_reset}                  Fasta file containing translated protein sequences from the GFF for running evaluation. If not specified, the workflow will automatically extract it from the `genome` and `gff`
    ${c_blue}--ref_protein${c_reset}              Fasta file containing the reference protein sequences to be used for evaluation. Ideally this should come from the same species and/or closely related species. If not provided, the Uniprot SwissProt data will be downloaded and used.
    ${c_blue}--lineage${c_reset}                  Lineage information providing for BUSCO, if not provided, the `--auto-lineage` option will be used instead. Example: eudicots_odb10
    ${c_blue}--oma_database${c_reset}             Pathway to the OMA database, if not specified, the workflow will download it automatically.
    ${c_blue}--run_blast${c_reset}                If specify, will use `blast` for running best reciprocal hits instead of DIAMOND. [default: false]

    ${c_yellow}Output option:${c_reset}
    ${c_blue}--outdir${c_reset}                   Output directory. [default: $params.outdir]
    ${c_blue}--tracedir${c_reset}                 Pipeline information. [default: $params.tracedir]
    ${c_blue}--publish_dir_mode${c_reset}         Option for nextflow to move data to the output directory. [default: $params.publish_dir_mode]
    ${c_blue}--tmpdir${c_reset}                   Database directory. [default: $params.tmpdir]

    ${c_yellow}Skipping options:${c_reset}
    ${c_dim}--skip_rename${c_reset}               Skip renaming genome fasta file by funannotate sort. [default: $params.skip_rename]

    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
    
    ${c_blue}Engines${c_reset} (choose one):
      docker
      singularity
      apptainer
    
    Per default: -profile slurm,singularity is executed. 
    ${c_reset}
    """.stripIndent()
}