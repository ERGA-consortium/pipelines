# Assembly Evaluation for ERGA-BGE Reports

The workflow requires the following:
- Species Taxonomy ID number
- NCBI Genome assembly accession code
- BUSCO Lineage
- WGS reads the accession code
- NCBI HiC reads accession code

It will obtain the data and process it to generate genome profiling (genomescope, smudgeplot is optional), assembly stats (gfastats), merqury stats (QV, completeness), BUSCO, snailplot, contamination blobplot, and HiC heatmap.

## One Assmebly, Illumina WGS reads + HiC reads 
![asm_rev_wf_illu](pics/asm_rev_1.png)

## One Assmebly, HiFi WGS reads + HiC reads
![asm_rev_wf_hifi](pics/asm_rev_2.png)
