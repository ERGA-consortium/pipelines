#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASSESS } from './workflows/assessment.nf'

workflow {
    println "START THE ANALYSIS\n"
    ASSESS()
}