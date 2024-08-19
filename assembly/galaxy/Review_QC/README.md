# Assembly Evaluation for ERGA-BGE Reports

The workflow requires the following:
- Species Taxonomy ID number
- NCBI Genome assembly accession code
- BUSCO Lineage
- WGS accurate reads accession code
- NCBI HiC reads accession code

The workflow will get the data and process it to generate genome profiling (genomescope, smudgeplot -optional-), assembly stats (gfastats), merqury stats (QV, completeness), BUSCO, snailplot, contamination blobplot, and HiC heatmap.

## One Assmebly, Illumina WGS reads + HiC reads
Use this workflow for ONT-based assemblies where the WGS accurate reads are Illumina PE
![asm_rev_wf_illu](pics/asm_rev_1.png)

## One Assmebly, HiFi WGS reads + HiC reads
Use this workflow for HiFi-based assemblies where the WGS accurate reads are PacBio HiFi
![asm_rev_wf_hifi](pics/asm_rev_2.png)
