# HiFi-Hic Assembly Pipeline


This pipeline serves as a collection of tools designed to facilitate genome assembly from HiFi PacBio data. This pipeline offers functionality for both local and server (slurm) deployment and it makes use of a series of conda environments. 
The workflow is structured into distinct units:
1. Preprocessing
- **status** Work in Progress
- **tools** cutadapt, HiFiadapterFilt, multiQC, fastQC
- **description** This unit involves trimming and Quality Control (fastQC and multiQC) of HiFi and HiC reads. Note that this phase is yet to be implemented.
2. Genome Profiling. 
- **tools** meryl, smudgeplot, genomescope
- **description**  Analyzes HiFi reads to extract critical information such as the estimated genome size and maximum depth coverage. This information is channeled to the next steps
3. Contigging. 
- **tools** Hifiasm, Canu, Flye
- **description**  This is the central part of the pipeline, where users can choose from multiple tools like Hifiasm, Canu and Flye (IMPORTANT the last two should still be tested, try useGrid=False, due to irregular behavior observed when running on the cluster). Hifiasm offers phasing mode with HiC data and various levels of purging (detailed in the Hifiasm parameter section).
4. Evaluation. 
- **tools** gfastats, busco, merqury
- **description**  Several tools are used to evaluate the assembly quality. Note: The integration of Busco is currently in progress.
5. Scaffolding 
- **status** Work in Progress
- **tools** YAHS
- **description** this unit will incorporate the tool YAHS when HiC data is provided. The reason is that we first need to decide how to choose (automatically? manually?) the purge level out of the three results of Hifiasm. 

![Hifi-HiC pipeline](https://github.com/valegale/pipelines/assets/74873652/9e1c3012-b8ee-4169-a699-0d39f66e7479)
