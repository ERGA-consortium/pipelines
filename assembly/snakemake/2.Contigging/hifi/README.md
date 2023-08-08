# HiFi-based contig assemblies

Currently performs contig assemblies from hifi reads based on 3 assemblers, hifiasm, canu and flye

Install the hifiasm and flye conda environments with the following commands:

```
conda create -n hifiasm -c bioconda -c conda-forge hifiasm

conda create -n canu -c bioconda -c conda-forge canu=2.2

conda create -n flye -c bioconda -c conda-forge flye
```

or:

```
conda create -n hifiasm -c bioconda -c conda-forge --file hifiasm.txt

conda create -n canu -c bioconda -c conda-forge --file canu.txt

conda create -n flye -c bioconda -c conda-forge --file flye.txt
```
