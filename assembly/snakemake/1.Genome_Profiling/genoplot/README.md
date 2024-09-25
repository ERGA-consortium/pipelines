# Genoplots

A snakemake pipeline which takes high-accuracy whole-genome reads as input (e.g. PacBio HiFI, Illumina WGS, 10x linked-reads), creates a kmer database and runs Genomescope and Smudgeplot to get estimates on genome size, heterozygosity, ploidy and read coverage

To run, you will need to create the following conda environments:

`conda create -n genoplots -c conda-forge -c bioconda meryl==1.3 genomescope2 smudgeplot`

`conda create -n kmc_smudge_pairs -c conda-forge kmc_smudge_pairs`

or:

`conda create -n genoplots -c conda-forge -c bioconda --file requirements_genoplots.txt`

`conda create -n kmc_smudge_pairs -c conda-forge kmc_smudge_pairs`

And edit the file names in the config.yaml file

You need to choose a kmer size to run all tools (default 21)

For genomescope2, you need to set the ploidy (default 2)

For smudgeplot, minimum and maximum thresholds for the kmer counts will be estimated by the kmer histogram as:

`conda activate genoplots && smudgeplot.py cutoff prefix.hist L`

to get an estimate of the lower bound. And

`conda activate genoplots && smudgeplot.py cutoff prefix.hist U`

to get an estimate of the upper bound.

