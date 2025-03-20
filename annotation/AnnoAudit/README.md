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
git clone https://github.com/ERGA-consortium/pipelines
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
      --taxon_id 9606 [Optional] \
      --ncbi_query_email xxxx \
      --rm -resume
```

### Other parameters for running the analysis

```
Input parameter:
--genome                  Draft genome fasta file contain the assembled contigs/scaffolds.
--gff                     Annotation file that need to be evaluated.
--genome_bam              BAM file contain the mapped information from the RNASeq to the genome FASTA.
--rnaseq                  A metadata CSV file following the pattern: sample_id,R1_path,R2_path,read_type. Required if `genome_bam` is not provided.
--taxon_id                Taxon ID for identifying BUSCO lineage and download protein data from NCBI if needed.
--ncbi_query_email        Email for querying protein from NCBI database.

Optional input:
--protein                  Fasta file containing translated protein sequences from the GFF for running evaluation. If not specified, the workflow will automatically extract it from the
 `genome` and `gff`.
--ref_protein              Fasta file containing the reference protein sequences to be used for evaluation. Ideally this should come from the same species and/or closely related specie
s. If not provided, the workflow will download the proteome from NCBI or using Uniprot SwissProt database.
--lineage                  Lineage information providing for BUSCO, if not provided, the workflow will automatically search for the closest lineage. Example: eudicots_odb10.
--genetic_code             Genetic code for translation of protein.
--stranding                Strandness of the RNASeq reads used for extraction of junction position using `regtools`.

Database input:
--odb_version              odb version to choose to run BUSCO, option: odb12, odb10. [default: odb12]
--busco_database           Pathway to the BUSCO databse store locally. [default: null]
--oma_database             Pathway to the OMA database, if not specified, the workflow will download it automatically. [default: null]
--ref_protein              Pathway to the reference proteome for comparison. [default: null]
--ncbi_query_count         Number of protein to extract from the NCBI database. [default: 100000]
--ncbi_query_batch         Number of protein to query for each batch. [default: 1000]

Output option:
--pdf                      Output PDF name. [default: AnnoAudit_Report.pdf]
--outdir                   Output directory. [default: /env/export/bigtmp2/pdoan/evaluate_pipeline]
--tracedir                 Pipeline information. [default: /env/export/bigtmp2/pdoan/evaluate_pipeline/pipeline_info]
--publish_dir_mode         Option for nextflow to move data to the output directory. [default: copy]
--tmpdir                   Database directory. [default: /env/export/bigtmp2/pdoan/evaluate_pipeline/tmpdir]

Conditioning options:
--rnaseq_single             If specify, will run `featureCounts` in single read mode, this is necessary if the mapped RNASeq is single-ended. [default: false]
--run_blast                 If specify, will use `blast` for running best reciprocal hits instead of DIAMOND. [default: false]
--query_ncbi_prot           If specify, will download the reference proteome from NCBI, other wise, will use the provided proteom or Uniprot SwissProt. [default: true]
--cds_only                  If specify, only extracting information from the GFF file using the CDS line. [default: "False"]

--help                   Print help message.

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

Below is the sample output of this workflow. The example PDF output is located in `assest` folder.

```
|General Statistics                 | Value           |
-------------------------------------------------------
|num_genes                          | 36391           |
|num_genes_without_introns          | 12968 (35.64%)  |
|mean_gene_length                   | 2359.57         |
|median_gene_length                 | 1562.0          |
|num_exons                          | 149725          |
|mean_exons_per_gene                | 4.11            |
|median_exons_per_gene              | 2.0             |
|num_exon_3n                        | 76783 (51.28%)  |
|num_exon_3n1                       | 36932 (24.67%)  |
|num_exon_3n2                       | 36010 (24.05%)  |
|mean_cds_length                    | 1091.4          |
|median_cds_length                  | 873.0           |
|total_cds_length                   | 39717145        |
|percentage_cds_coverage            | 10.64%          |
|num_introns                        | 113334          |
|mean_intron_length                 | 407.2           |
|median_intron_length               | 149.0           |
|short_intron_<120_3n0_without_stop | 4324 (3.82)%    |
|long_intron_>120_3n0_without_stop  | 1185 (1.05)%    |
|short_intron_<120_3n1_without_stop | 4205 (3.71)%    |
|long_intron_>120_3n1_without_stop  | 1291 (1.14)%    |
|short_intron_<120_3n2_without_stop | 4319 (3.81)%    |
|long_intron_>120_3n2_without_stop  | 1249 (1.10)%    |
|short_intron_<120_3n0_with_stop    | 12073 (10.65)%  |
|long_intron_>120_3n0_with_stop     | 20332 (17.94)%  |
|short_intron_<120_3n1_with_stop    | 11652 (10.28)%  |
|long_intron_>120_3n1_with_stop     | 20486 (18.08)%  |
|short_intron_<120_3n2_with_stop    | 11733 (10.35)%  |
|long_intron_>120_3n2_with_stop     | 20485 (18.07)%  |

|BUSCO                              | Value           |
-------------------------------------------------------
|lineage_dataset                    | poales_odb10    |
|complete                           | 97.6%           |
|single_copy                        | 95.8%           |
|multi_copy                         | 1.8%            |
|fragmented                         | 0.2%            |
|missing                            | 2.2%            |
|num_markers                        | 4896            |
|domain                             | eukaryota       |


|OMARK                              | Value           |
-------------------------------------------------------
|OMA_clade                          | Oryza           |
|num_conserved_hogs                 | 15087           |
|single                             | 13316 (88.26%)  |
|duplicated                         | 1353 (8.97%)    |
|duplicated_unexpected              | 1101 (7.30%)    |
|duplicated_expected                | 252 (1.67%)     |
|missing                            | 418 (2.77%)     |
|num_proteins_in_proteome           | 36387           |
|total_consistent                   | 30365 (83.45%)  |
|consistent_partial_hits            | 1803 (4.96%)    |
|consistent_fragmented              | 1625 (4.47%)    |
|total_inconsistent                 | 2283 (6.27%)    |
|inconsistent_partial_hits          | 517 (1.42%)     |
|inconsistent_fragmented            | 1444 (3.97%)    |
|total_contaminants                 | 0 (0.00%)       |
|contaminants_partial_hits          | 0 (0.00%)       |
|contaminants_fragmented            | 0 (0.00%)       |
|total_unknown                      | 3739 (10.28%)   |

|PSAURON                            | Value           |
-------------------------------------------------------
|psauron_score                      | 83.8            |
|true_count                         | 30494           |
|false_count                        | 5893            |
|median_score                       | 0.98278         |
|max_score                          | 1.0             |
|min_score                          | 0.00022         |


|Best Reciprocal Hits               | Value           |
-------------------------------------------------------
|num_best_reciprocal_hits           | 29185           |
|num_splitting_genes_08             | 932 (3.19%)     |
|num_splitting_genes_05             | 0 (0.0%)        |
|num_fusion_genes_12                | 437 (1.5%)      |
|num_fusion_genes_15                | 482 (1.65%)     |
|KL_divergence_normalized           | 0.0105          |
|JS_divergence_normalized           | 0.0023          |
|Wasserstein_distance               | 2.480915        |


|RNASeq                             | Value           |
-------------------------------------------------------
|mapping_rate                       | 96.27%          |
|primary_mapping_rate               | 95.83%          |
|properly_paired                    | 92.47%          |
|num_gene_unsupported               | 9445 (25.95%)   |
|num_exon_unsupported               | 20232 (13.51%)  |
|num_intron_supported               | 107202          |
|num_intron_supported_canonical     | 107131 (99.93%) |
|num_intron_supported_non_canonical | 71 (0.07%)      |
```

## Performance of the workflow on assessing annotation

To be added

## Future work

- Adding other plots for easier evaluation
- Perform comparative performance with different genomes