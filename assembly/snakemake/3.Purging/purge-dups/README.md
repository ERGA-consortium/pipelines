## Purge dups

This snakemake pipeline is designed to be run using as input a contig-level genome and pacbio reads. To create the necessary conda environment run:

conda create --name purge-dups --file requirements.txt

Run in terminal with command:

snakemake --use-conda --cores N

Or configure the cluster.json and run using the `run_cluster` command
