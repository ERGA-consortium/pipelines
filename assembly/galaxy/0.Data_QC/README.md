## Reads preprocessing & QC
Galaxy Workflows for checking quality and trimming reads from Illumina (WGS, HiC), PacBio and ONT.

Load the respective .ga file in Galaxy to run the workflow.


#### Illumina Reads
![QCillu2305](pics/QCillu2305.png)
The workflow takes a [paired-reads collection](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html) and runs FastQC and Seqkit, trims with Fastp, and creates a MultiQC report. The main outputs are a paired collection of trimmed reads, a report with raw and trimmed reads stats, and a table with raw reads stats.

#### PacBio Reads


#### ONT Reads
![QCont2305](pics/QCont2305.png)
The workflow takes ONT reads collection and runs SeqKit, filters by quality and plots the result. The main outputs are a collection of filtered reads, a table with raw reads stats and a plot of filtered reads.


