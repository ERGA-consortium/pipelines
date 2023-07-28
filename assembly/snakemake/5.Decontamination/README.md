# Decontamination

This pipeline runs Blobtools to identify potential contaminants in the assembly. By default any non-eukaryotic sequences are removed, but the blobtools outputs allow for further manual inspection in case of other contaminations in the assembly.

Create the conda environment using either of the following commands:

`conda create -n decontamination -c bioconda -c conda-forge blast blobtools minimap2 samtools bbmap`

or:

`conda create -n decontamination -c bioconda -c conda-forge --file decontamination.txt`

You will need to download the nt database and store this locally. The directory containing the nt database should be listed in the `config.yaml` files under `nt_dbDir`

To download/update the database run:

```cd /db/nt/

conda activate decontamination

update_blastdb.pl --source aws -decompress --passive --num_threads {threads} blastdb nt

conda deactivate```

This will take some time and may restart if the download fails at any point.

Create a file of files which lists the long-read sequences used to assemble your genome. These will be used to determine the read-coverage of the sequences in your assembly. For example:

`ls /data/reads/pacbio/*fasta > lr.fofn`

or:

`ls /data/reads/ont/*fastq.gz > lr.fofn'

