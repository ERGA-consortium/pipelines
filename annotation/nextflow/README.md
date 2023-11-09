# ANNOTATO - Annotation workflow To Annotate Them Oll

## Prerequisites

The following programs are required to run the workflow and the listed version were tested. 

`nextflow v23.04.0 or higher`

`conda` or `singularity`

`docker` (have not been tested)

## Workflow

The pipeline is based on `Funannotate` and/or `BRAKER` and was tested on these following datasets: 
- Drosophila melanogaster: [https://doi.org/10.5281/zenodo.8013373](https://doi.org/10.5281/zenodo.8013373)
- Pocillopora meandrina

### Input data

- Reference genome `genome.[.fna, .fa, .fasta][.gz]`
- RNAseq data listed in a metadata csv file. Input type can be mixed between long and short reads, with the option of single-end read. The input file should follow the format below:

```
sample_id,R1_path,R2_path,read_type
SAM1,/path/to/R1,,long             # For long reads
SAM2,/path/to/R1,/path/to/R2,short # For PE reads
SAM3,/path/to/R1,,short            # For SE reads
```

- Protein sequence data in fasta format, could be gzip or not

### Pipeline steps

![Pipeline](./assets/images/annotato-workflow.drawio.svg)

The main pipeline is divided into five different subworkflows.
- `Preprocess RNA` is where the input RNASeq data are QC and trimmed.
- `Process RNA Minimap` is triggered when long reads FastQ are in the input CSV file.
- `Process RNA STAR` will run when short reads FastQ are in the input CSV.
- `Genome Masking` runs by default if not skipped. It assumes the input genome fasta is not masked and will run Denovo repeat masking with RepeatModeler and RepeatMasker.
- `Filter Repeat` whenever there is a Denovo masking step, this sub-workflow will be triggered to remove the repeat sequences that appeared in the Uniprot Swissprot protein data. 

### Output data

- MultiQC report for the RNASeq data, before and after trimming, mapping rate of short reads, and the BUSCO results of predicted genes.
- RepeatMasker report containing quantity of masked sequence and distribution among TE families
- Protein-coding gene annotation file in gff3 format
- BUSCO summary of annotated sequences

## Running Annotato

### Before running the pipeline

One thing with Nextflow is that it is running off a Java Virtual Machine (JVM), and it will try to use all available memory for Nextflow even though it is unnecessary (for workflow management and job control). This will cause much trouble if you run a job on an HPC cluster. Thus, to minimize the effect of it, we need to limit the maximum memory the JVM can use.

```
export NFX_OPTS="-Xms=512m -Xmx=3g"
```

`-Xms` is the lower limit, which is set as 512 MB.
`-Xmx` is the upper limit, which in this case is set as 3 GB.
Please modify this according to your situation.

### Without RNASeq and protein data

Perform the analysis with only the draft genome and busco database.

```
nextflow run main.nf --genome /path/to/genome.fasta --species "Abc def" --buscodb 'metazoa' 
```

The workflow will run Denovo repeat masking on the draft genome, then softmask the repeat region and use the genome to run `funannotate`. Add `--run_braker` to run the genome prediction using `BRAKER` instead.

### Running Annotato with RNASeq data

When you want to let the workflow run the mapping by itself, uses `input.csv` as input with the link to all `FASTQ` file.

```
nextflow run main.nf --genome /path/to/genome.fasta[.gz] --rnaseq /path/to/input.csv --species "Abc def" --buscodb 'metazoa' 
```

Based on the content of the `input.csv` file to trigger different RNASeq processing workflows. The output `bam` file will then be used for genome prediction.

When reads are mapped to the reference genome, the aligned `bam` file can be used as input to the pipeline instead of the raw `FASTQ`

```
nextflow run main.nf --genome /path/to/genome.fasta[.gz] --short_rna_bam /path/to/shortreads.bam [--long_rna_bam /path/to/longreads.bam] --species "Abc def" --buscodb 'metazoa' 
```

**ATTENTION**: One major drawback of the current workflow is that the input genome will be sorted and renamed by the `funannotate sort` function. This is because `AUGUSTUS` and `Funannotate` won't work normally when the header of the input genome is too long and contains weird characters. Therefore, if you want to provide a `bam` file as input instead of the raw `FASTQ`, please run `funannotate sort` on the genome fasta first and then use it as the reference for running alignment. Or in case your genome headers are already shorter than 16 character, please add `--skip_rename` when running the pipeline.

### Running Annotato with protein data

```
nextflow run main.nf --genome /path/to/genome.fasta[.gz] --protein /path/to/protein.fasta[.gz] --species "Abc def" --buscodb 'metazoa' 
```

When only protein data is provided, the workflow will run denovo masking then repeat filter with the additional protein data. The masked genome and protein fasta will then be used for gene prediction.

### Running Annotato with both protein and RNASeq data

The full pipeline is triggered when both RNASeq data and protein fasta is provided.

```
nextflow run main.nf --genome /path/to/genome.fasta[.gz] --protein /path/to/protein.fasta[.gz] --rnaseq /path/to/input.csv --species "Abc def" --buscodb 'metazoa' 
```

### Other parameters for running the analysis

```
Compulsory input:
--genome                       Draft genome fasta file contain the assembled contigs/scaffolds
--species                      Species name for the annotation pipeline, e.g. "Drosophila melanogaster"

Optional input:
--protein                      Fasta file containing known protein sequences used as an additional information for gene prediction pipeline.
                               Ideally this should come from the same species and/or closely related species. [default: null]
--rnaseq                       A CSV file following the pattern: sample_id,R1_path,R2_path,read_type.
                               This could be generated using gen_input.py. Run `python gen_input.py --help` for more information. 
                               [default: null]
--long_rna_bam                 A BAM file for the alignment of long reads (if any) to the draft genome. Noted that the header of the draft
                               genome need to be renamed first before alignment otherwise it will causes trouble for AUGUSTUS and funannotate. 
                               [default: null]
--short_rna_bam                A BAM file for the alignment of short reads (if any) to the draft genome. Noted that the header of the draft 
                               genome need to be renamed first before alignment otherwise it will causes trouble for AUGUSTUS and funannotate. 
                               [default: null]
--knownrepeat                  Fasta file containing known repeat sequences of the species, this will be used directly for masking 
                               (if --skip_denovo_masking) or in combination with the denovo masking. [default: null]

Output option:
--outdir                       Output directory. 
--tracedir                     Pipeline information. 
--publish_dir_mode             Option for nextflow to move data to the output directory. [default: copy]
--tmpdir                       Database directory. 

Funannotate params:
--run_funannotate              Whether to use funannotate for gene prediction. [default: true]
--organism                     Fungal-specific option. Should be change to "fungus" if the annotated organism is fungal. [default: other]
--ploidy                       Set the ploidy for gene prediction, in case of haploid, a cleaning step will be performed by funannotate to remove
                               duplicated contigs/scaffold. [default: 2]
--buscodb                      BUSCO database used for AUGUSTUS training and evaluation. [default: eukaryota]
--buscoseed                    AUGUSTUS pre-trained species to start BUSCO. Will be override if rnaseq data is provided. [default: null]

Skipping options:
--skip_rename                  Skip renaming genome fasta file by funannotate sort. 
--skip_all_masking             Skip all masking processes, please be sure that your --genome input is soft-masked before triggering this 
                               parameter. [default: false]
--skip_denovo_masking          Skip denovo masking using RepeatModeler, this option can only be run when --knownrepeat fasta is provided. 
                               [default: false]
--skip_functional_annotation   Skip functional annotation step. [default: false]
--skip_read_preprocessing      Skip RNASeq preprocessing step. [default: false]

Execution/Engine profiles:
The pipeline supports profiles to run via different Executers and Engines e.g.: -profile local,conda

Executer (choose one):
  local
  slurm

Engines (choose one):
  conda
  mamba
  docker
  singularity

Per default: -profile slurm,singularity is executed.
```

## Future work
- Adding functional annotation with `Interproscan` and `eggnog`