# AnnoAudit - Annotation Auditor

AnnoAudit is a robust Nextflow pipeline designed to evaluate the quality of genomic annotations through a multifaceted approach.

## Overview of the workflow

The workflow assess the annotation quality based on different criteria:
- Protein evidence support
- RNASeq evidence support
- Statistics of the predictions (i.e., gene length, exon number, etc.)
- Ortholog analysis (BUSCO, OMArk)

### Input data

- Reference genome `genome.[.fna, .fa, .fasta]`
- Annotation output `annotation.gff`
- RNAseq data listed in a metadata csv file. Input type can be mixed between long and short reads, with the option of single-end read. The input file should follow the format below:

```
sample_id,R1_path,R2_path,read_type
SAM1,/path/to/R1,,long             # For long reads
SAM2,/path/to/R1,/path/to/R2,short # For PE reads
SAM3,/path/to/R1,,short            # For SE reads
```

- Protein reference data in `fasta` format for evaluation, if not given, then the `Uniprot-SwissProt` will be downloaded and used.

### Pipeline steps

![Pipeline](./assets/images/annoaudit-workflow.svg)

The main pipeline is divided into five different subworkflows.
- `General statistics`: Calculate the statistics obtained from the GFF file.
- `RNASeq analysis`: Map the RNASeq data to the genome (or with provided mapping bam file) to generate exon, intron, transcript coverage.
- `Ortholog analysis`: Compare the predicted proteome to known database using BUSCO and OMArk (OMA database).
- `Protein analysis`: Blast the predicted proteome to a known database (could be of relative species) to obtain best reciprocal hits (BRH), then generate statistics based on the BRH results.

### Output data

- Output text file contain the statistic calculated from the input `GFF` file: 
  - General statistics
  - BUSCO
  - OMArk
  - Best reciprocal hits
  - RNASeq analysis

## Prerequisites

The following programs are required to run the workflow and the listed version were tested. 

`nextflow v23.04.0 or higher`

`singularity`

`docker` (have not been tested but in theory should work fine)

## Installation

Simply get the code from github or workflowhub and directly use it for the analysis with `nextflow`.

```
git clone https://github.com/phuongdoand/ERGA-pipelines/tree/main/annotation/annotation_assessment
```

## Running AnnoAudit

### Before running the pipeline (IMPORTANT)

One thing with Nextflow is that it is running off a Java Virtual Machine (JVM), and it will try to use all available memory for Nextflow even though it is unnecessary (for workflow management and job control). This will cause much trouble if you run a job on an HPC cluster. Thus, to minimize the effect of it, we need to limit the maximum memory the JVM can use.

```
export NFX_OPTS="-Xms=512m -Xmx=3g"
```

`-Xms` is the lower limit, which is set as 512 MB.
`-Xmx` is the upper limit, which in this case is set as 3 GB.
Please modify this according to your situation.

### How to run the code

```
nextflow run main.nf --genome genome.fasta \
      --gff annotation.gff3 \
      --rnaseq metadata.csv [--genome_bam path/to/the/mapped/bam]\
      --outdir OUTDIR_NAME \
      --rm -resume
```

### Other parameters for running the analysis

```
Compulsory input:
--genome                       Draft genome fasta file contain the assembled contigs/scaffolds
--gff                          Annotation file that need to be evaluated
--genome_bam                   BAM file contain the mapped information from the RNASeq to the genome FASTA, will be utilize instead of 
                               `--rnaseq` once given.
--rnaseq                       A CSV file following the pattern: sample_id,R1_path,R2_path,read_type.
--genetic_code                 Genetic code for translate the intron sequences. [default: 1]

Optional input:
--protein                      Fasta file containing translated protein sequences from the GFF for running evaluation.
                               If not specified, the workflow will automatically extract it from the `genome` and `gff`
--ref_protein                  Fasta file containing the reference protein sequences to be used for evaluation.
                               Ideally this should come from the same species and/or closely related species. If not provided
                               the workflow will download the protein sequences from the NCBI Entrez database (priority), if skipped as well, the Uniprot SwissProt data will be downloaded and used.
--lineage                      Lineage information providing for BUSCO, if not provided, the `--auto-lineage` option will be used
                               instead. Example: eudicots_odb10
--oma_database                 Pathway to the OMA database, if not specified, the workflow will download it automatically.

Query NCBI option:
--taxon_id                     Taxon ID required for searching the NCBI database, required if the `--query_ncbi_prot` is set (by default), 
                               and the `--ref_protein` is not provided.
--ncbi_query_email             Email address to use for NCBI Entrez, required if the `--query_ncbi_prot` is set (by default), 
                               and the `--ref_protein` is not provided.
--ncbi_query_count             Number of protein to query from the Entrez database. [default: 100000]
--ncbi_query_batch             Number of protein to search at a time. [default: 1000]

Running option:
--cds_only                     Whether or not to extract the exon data from the CDS line in the provided GFF. [default: False] [option: True, False]
--query_ncbi_prot              Query protein data from the NCBI Entrez database. [default: true]
--run_blast                    If specify, will use `blast` for running best reciprocal hits instead of DIAMOND. [default: false]
--skip_omark                   Skip OMArk analysis step. [default: false] (to be added in the future)

Output option:
--outdir                       Output directory. 
--tracedir                     Pipeline information. 
--publish_dir_mode             Option for nextflow to move data to the output directory. [default: copy]
--tmpdir                       Database directory. 

Execution/Engine profiles:
The pipeline supports profiles to run via different Executers and Engines e.g.: -profile local,conda

Executer (choose one):
  local
  slurm

Engines (choose one):
  docker
  singularity
  apptainer

Per default: -profile slurm,singularity is executed.
```

## Example output

Below is the sample output of this workflow

```
|General Statistics                 | Value                |
------------------------------------------------------------
|num_genes                          | 41048                |
|num_genes_without_introns          | 14365 (35.0%)        |
|mean_gene_length                   | 2444.24              |
|median_gene_length                 | 1613.0               |
|num_exons                          | 175037               |
|mean_exons_per_gene                | 4.26                 |
|median_exons_per_gene              | 2.0                  |
|num_exon_3n                        | 89250 (50.99%)       |
|num_exon_3n1                       | 43323 (24.75%)       |
|num_exon_3n2                       | 42464 (24.26%)       |
|mean_cds_length                    | 1109.73              |
|median_cds_length                  | 894.0                |
|total_cds_length                   | 45552283             |
|percentage_cds_coverage            | 12.2%                |
|num_introns                        | 133989               |
|mean_intron_length                 | 408.83               |
|median_intron_length               | 147.0                |
|short_intron_<120_3n0_without_stop | 5027 (3.75)%         |
|long_intron_>120_3n0_without_stop  | 1503 (1.12)%         |
|short_intron_<120_3n1_without_stop | 5027 (3.75)%         |
|long_intron_>120_3n1_without_stop  | 1515 (1.13)%         |
|short_intron_<120_3n2_without_stop | 5126 (3.83)%         |
|long_intron_>120_3n2_without_stop  | 1473 (1.10)%         |
|short_intron_<120_3n0_with_stop    | 13849 (10.34)%       |
|long_intron_>120_3n0_with_stop     | 24363 (18.18)%       |
|short_intron_<120_3n1_with_stop    | 13849 (10.34)%       |
|long_intron_>120_3n1_with_stop     | 24076 (17.97)%       |
|short_intron_<120_3n2_with_stop    | 14068 (10.50)%       |
|long_intron_>120_3n2_with_stop     | 24113 (18.00)%       |

|BUSCO                              | Value                |
------------------------------------------------------------
|lineage_dataset                    | eudicotyledons_odb10 |
|complete                           | 86.7%                |
|single_copy                        | 70.5%                |
|multi_copy                         | 16.2%                |
|fragmented                         | 0.7%                 |
|missing                            | 12.6%                |
|num_markers                        | 2326                 |
|domain                             | eukaryota            |

|OMARK                              | Value                |
------------------------------------------------------------
|OMA_clade                          | Oryza                |
|num_conserved_hogs                 | 15087                |
|single                             | 11686 (77.46%)       |
|duplicated                         | 2986 (19.79%)        |
|duplicated_unexpected              | 2749 (18.22%)        |
|duplicated_expected                | 237 (1.57%)          |
|missing                            | 415 (2.75%)          |
|num_proteins_in_proteome           | 41044                |
|total_consistent                   | 34598 (84.29%)       |
|consistent_partial_hits            | 1973 (4.81%)         |
|consistent_fragmented              | 1977 (4.82%)         |
|total_inconsistent                 | 2456 (5.98%)         |
|inconsistent_partial_hits          | 567 (1.38%)          |
|inconsistent_fragmented            | 1541 (3.75%)         |
|total_contaminants                 | 0 (0.00%)            |
|contaminants_partial_hits          | 0 (0.00%)            |
|contaminants_fragmented            | 0 (0.00%)            |
|total_unknown                      | 3990 (9.72%)         |

|Best Reciprocal Hits               | Value                |
------------------------------------------------------------
|num_best_reciprocal_hits           | 29076                |
|num_splitting_genes_08             | 1578 (5.43%)         |
|num_splitting_genes_05             | 0 (0.0%)             |
|num_fusion_genes_12                | 331 (1.14%)          |
|num_fusion_genes_15                | 342 (1.18%)          |
|KL_divergence_normed               | 0.0186               |
|JS_divergence_normed               | 0.0039               |
|Wasserstein_distance               | 8.469494             |

|RNASeq                             | Value                |
------------------------------------------------------------
|mapping_rate                       | 96.27%               |
|primary_mapping_rate               | 95.83%               |
|properly_paired                    | 92.47%               |
|num_gene_unsupported               | 15472 (37.69%)       |
|num_exon_unsupported               | 23003 (13.14%)       |
```

## Performance of the workflow on assessing annotation

To be added

## Future work

- Adding more criteria for the RNASeq analysis
- Perform comparative performance with different genomes
