# YAHS scaffolding

A snakemake pipeline to perform scaffolding using Chromatin Conformation Capture (HiC) reads with the software YAHS.

To run, you will need to create the following conda environment:

`conda create -n yahs -c conda-forge -c bioconda yahs pairtools bwa samtools numpy==1.23`

or:

`conda create -n yahs -c conda-forge -c bioconda --file requirements_yahs.txt`


