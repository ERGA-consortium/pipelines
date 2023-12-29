# HiFi-Hic Assembly Pipeline


This pipeline serves as a collection of tools designed to facilitate genome assembly from HiFi PacBio data. This pipeline offers functionality for both local and server (slurm) deployment and it makes use of a series of conda environments. 
The workflow is structured into distinct units:
1. Preprocessing
- **status** Work in Progress
- **tools** cutadapt, HiFiadapterFilt, multiQC, fastQC
- **description** This unit involves trimming and Quality Control (fastQC and multiQC) of HiFi and HiC reads. Note that this phase is yet to be implemented.
2. Genome Profiling 
- **tools** meryl, smudgeplot, genomescope
- **description**  Analyzes HiFi reads to extract critical information such as the estimated genome size and maximum depth coverage. This information is channeled to the next steps
3. Contigging
- **tools** Hifiasm, Canu, Flye
- **description**  This is the central part of the pipeline, where users can choose from multiple tools like Hifiasm, Canu and Flye (IMPORTANT the last two should still be tested, try useGrid=False, due to irregular behavior observed when running on the cluster). Hifiasm offers phasing mode with HiC data and various levels of purging (detailed in the Hifiasm parameter section).
4. Evaluation
- **tools** gfastats, busco, merqury
- **description**  Several tools are used to evaluate the assembly quality. Note: The integration of Busco is currently in progress.
5. Scaffolding 
- **status** Work in Progress
- **tools** YAHS
- **description** this unit will incorporate the tool YAHS when HiC data is provided. Still not implemented since I am not sure how to choose (automatically? manually?) the purge level out of the three results of Hifiasm with the "all" option for purgeL. 

![Hifi-HiC pipeline](https://github.com/valegale/pipelines/assets/74873652/9e1c3012-b8ee-4169-a699-0d39f66e7479)


## Running the pipeline


To execute the pipeline, ensure that [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) are correctly installed on your system.
The pipeline can be launched with the following commands:

#### Default Execution
```
snakemake --profile submit_config/default/
```

This command runs the pipeline with default parameters: 

```
cores: 32
dry-run: False
use-conda: True
latency-wait: 60
rerun-incomplete: True
printshellcmds: False
keep-going: False
```


#### HPC (slurm) submission

```
snakemake --profile submit_config/slurm/
```

Use this command to submit the pipeline via an HPC system. Resource specifications for each tool are defined in `submit_config/slurm/define_resources.yaml`.

### Config parameters

Paths to Input/Output Files and Folders
- **input_file_hifi.tsv**  Path(s) to the file(s) containing HiFi data paths. Each row corresponds to paths for the same species, within a single column named **hifi_data** (refer to the sample file inside the **test** folder).
- **input_file_hic.tsv** Path to the file listing two Hi-C file paths. This file should have columns named **HiC_R1** and **HiC_R2** (refer to the sample file inside the **test** folder).
- **results** Path to the results folder.

General parameters
- **prefix** Prefix for the species under analysis used in output files/folders.
- **ploidy** Ploidy of the genome. 
- **kmersize** K-mer size used for counting k-mers in the reads.
- **trimming** Boolean (True/False) to enable/disable trimming.

Contigging parameters
- **contig_assembler** Options: hifiasm, flye, canu
- **purgel**
Hifiasm purging parameter:
0 = No purging
1 = Purge contained haplotigs (selectivity threshold s=0.75)
2 = Purge all haplotigs (selectivity threshold s=0.75)
3 = Aggressive purging of all haplotigs (selectivity threshold s=0.55)
"all" = sequentially execute purging steps 1, 2, and 3 (except zero).
- **phasing_mode** Valid only when `contig_assembler = hifiasm` and Hi-C files are provided. If `purgel = "all"`, phasing mode is automatically enabled.
 
Evaluation parameters
- **busco_lineage** Busco lineage specification (busco is not yet implemented).

Resource configuration
- **resources** Path to define Slurm resources.
