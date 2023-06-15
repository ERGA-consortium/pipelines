### Reads preprocessing & QC

Galaxy Workflows for checking quality and trimming reads from Illumina (WGS, HiC), PacBio and ONT.

Load the respective .ga file in Galaxy to run the workflow.

#### Illumina Reads
![QCillu2305](pics/QCillu2305.png)
The workflow takes a [paired-reads collection](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html) and runs FastQC and Seqkit, trimms with Fastp, and create a MultiQC report. The main outputs are a paired-collection of trimmed reads, a report with raw and trimmed reads stats, and a table with raw reads stats.

#### PacBio Reads


#### ONT Reads


