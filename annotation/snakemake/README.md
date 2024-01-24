Protein-coding gene annotation workflow. designed, developed, and tested by Sagane Joye (ERGA-CH, UNIL)

https://github.com/sdind/genome_annotation_workflow

## Prerequisites

The following programs are required to run the workflow and the listed version were tested. It should be noted that older versions of snakemake are not compatible with newer versions of singularity as is noted here: [https://github.com/snakemake/snakemake/issues/2319](https://github.com/snakemake/snakemake/issues/2319).

`conda v 23.7.3`
`singularity v 3.7.3`
`snakemake v 7.32.3` 

## Workflow

The pipeline is based on braker3 and was tested on the following dataset from Drosophila melanogaster: [https://doi.org/10.5281/zenodo.8013373](https://doi.org/10.5281/zenodo.8013373)

### Input data

- Reference genome in fasta format

- RNAseq data in paired-end zipped fastq format

- uniprot fasta sequences in zipped fasta format

### Pipeline steps

- **Repeat Model and Mask** Run RepeatModeler using the genome as input, filter any repeats also annotated as protein sequences in the uniprot database and use this filtered libray to mask the genome with RepeatMasker

- **Map RNAseq data** Trim any remaining adapter sequences and map the trimmed reads to the input genome

- **Run gene prediction software** Use the mapped RNAseq reads and the uniprot sequences to create hints for gene prediction using Braker3 on the masked genome

- **Evaluate annotation** Run BUSCO to evaluate the completeness of the annotation produced

### Output data

- FastQC reports for input RNAseq data before and after adapter trimming

- RepeatMasker report containing quantity of masked sequence and distribution among TE families

- Protein-coding gene annotation file in gff3 format

- BUSCO summary of annotated sequences

## Setup

Your data should be placed in the `data` folder, with the reference genome in the folder `data/ref` and the transcript data in the foler `data/rnaseq`.

The config file requires the following to be given:

```
asm: 'absolute path to reference fasta'
snakemake_dir_path: 'path to snakemake working directory'
name: 'name for project, e.g. mHomSap1'
RNA_dir: 'absolute path to rnaseq directory'
busco_phylum: 'busco database to use for evaluation e.g. mammalia_odb10'
```
