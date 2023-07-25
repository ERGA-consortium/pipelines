from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module assembly_workflow:
  snakefile: "../modules/assemble_genomes.rules.smk"

working_dir = config["Outputs"]["base_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

assembled = []

##0. Define path for files and variables

flye_dir = config["Outputs"]["flye_dir"]
nextdenovo_dir = config["Outputs"]["nextdenovo_dir"]
flye_assembly = config["Outputs"]["flye_out"]
nextdenovo_assembly = config["Outputs"]["nextdenovo_out"]


ONT_filtered = config["Inputs"]["ONT_filtered"] 

##1. Run assemblers
if config["Parameters"]["run_flye"] == True:
  assembled.append(flye_assembly)
  use rule flye from assembly_workflow with:
    input:
      reads = ONT_filtered,
    output:
      assembly = flye_assembly,
      gfa = report(flye_dir + "assembly_graph.gfa.png",
            caption="../report/flye.rst",
            category = "Evaluate assemblies",
            subcategory = flye_dir)
    params:
      outdir = flye_dir,
      readtype = config["Parameters"]["lr_type"],
      pol_iterations = config["Flye"]["Flye polishing iterations"],
      other_flye_opts = config["Flye"]["options"]
    log:
      flye_dir + "logs/" + str(date) + ".j%j.flye.out",
      flye_dir + "logs/" + str(date) + ".j%j.flye.err"
    benchmark:
      flye_dir + "logs/" + str(date) + ".flye.benchmark.txt",
    conda:
      '../envs/flye2.9.1.yaml'
    threads: config["Flye"]["Flye cores"]

if config["Parameters"]["run_nextdenovo"] == True:
 assembled.append(nextdenovo_assembly)
 use rule nextdenovo from assembly_workflow with:
  input:
    reads = ONT_filtered
  output:
    assembly = nextdenovo_assembly
  params:
    outdir = nextdenovo_dir,
    config = config["Parameters"]["ndconfFile"]
  log:
    nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.out",
    nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.err"
  benchmark:
    nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.benchmark.txt"
  envmodules:
    "NextDenovo/2.5.0"
  threads: config["Nextdenovo"]["Nextdenovo cores"]