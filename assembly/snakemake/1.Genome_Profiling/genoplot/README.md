# Genoplots

A snakemake pipeline which takes high-accuracy whole-genome reads as input (e.g. PacBio HiFI, Illumina WGS, 10x linked-reads), creates a kmer database and runs Genomescope and Smudgeplot to get estimates on genome size, heterozygosity, ploidy and read coverage

To run, you will need to create the following conda environment:

`conda create -n genoplots -c conda-forge -c bioconda meryl==1.3 genomescope2 smudgeplot`

or:

`conda create -n genoplots -c conda-forge -c bioconda --file requirements_genoplots.txt`

And edit the file names in the config.yaml file

You need to choose a kmer size to run all tools (default 21)

For genomescope2, you need to set the ploidy (default 2)

For smudgeplot, you need to set minimum and maximum thresholds for the kmer counts. By default these are set to a minimum of 10 and maximum of 1500. If coverage is very high or low, these should be adjusted. A good estimate can be calculated by running 

`conda activate genoplots && smudgeplot.py cutoff prefix.hist L`

to get an estimate of the lower bound. And

`conda activate genoplots && smudgeplot.py cutoff prefix.hist U`

to get an estimate of the upper bound. This needs to be run after meryl has created the .hist file, which is the input for genomescope2.

