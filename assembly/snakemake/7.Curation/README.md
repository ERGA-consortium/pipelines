# HiC contact map generation

Snakemake pipeline for the generation of `.pretext` and `.mcool` files for visualisation of HiC contact maps with the softwares PretextView and HiGlass, respectively.

## Prerequisites

This pipeine has been tested using `Snakemake v7.32.4` and requires conda for installation of required tools. To run the pipline use the command:

`snakemake --use-conda`

There are provided a set of configuration and running scripts for exectution on a slurm queueing system. After configuring the `cluster.json` file run:

`./run_cluster`

## Before starting

You need to create a temporary folder and specify the path in the `config.yaml` file. This should be able to hold the temporary files created when sorting the `.pairsam` file (100s of GB or even many TBs)

The path to the genome assemly must be given in the `config.yaml`.

The HiC reads should be paired and named as follows: `Library_1.fastq.gz Library_2.fastq.gz`. The pipeline can accept any number of paired HiC read files, but the naming must be consistent. The folder containing these files must be provided in the `config.yaml`.
