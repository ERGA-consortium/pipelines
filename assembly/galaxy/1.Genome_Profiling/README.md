## Kmer-based Genome Profiling
Galaxy Workflows for checking properties of the genome (size, heterozygosity, ploidy) based on reads (Illumina, HiFi) k-mer profile.

Load the respective .ga file in Galaxy to run the workflow.

### Illumina Reads
The workflow takes a trimmed paired-reads collection, runs Meryl to create a K-mer database, Genomescope2 to estimate genome properties and Smudgeplot to estimate ploidy. The main results are K-mer database and genome profiling plots and tables.
Default K-mer length and ploidy for Genomescope are 21 and 2, respectively.
![ProfIllu](pics/ProfIllu2305.png)

### HiFi Reads
