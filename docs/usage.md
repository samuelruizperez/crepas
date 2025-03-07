# grothlab/glseq: Usage

> [!IMPORTANT]
> Please read this documentation on the grothlab/glseq repository: [https://github.com/grothlab/glseq/blob/dev/docs/usage.md](https://github.com/grothlab/glseq/blob/dev/docs/usage.md)

## Table of Contents

1. [Samplesheet input](#samplesheet-input)
2. [Reference genome files](#reference-genome-files)
3. [Running the pipeline](#running-the-pipeline)
    1. [Parameters](#parameters)
    2. [Updating the pipeline](#updating-the-pipeline)
    3. [Reproducibility](#reproducibility)
    4. [Core Nextflow arguments](#core-nextflow-arguments)
    5. [Custom configuration](#custom-configuration)
    6. [Running in the background](#running-in-the-background)
    7. [Nextflow memory requirements](#nextflow-memory-requirements)
4. [**Quick start guide for DAN System users**](#quick-start-guide-for-dan-system-users)

---

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the `--input` parameter to specify its location. It has to be a comma-separated (`.csv`) file with with 11 columns and a header row as explained below.

```bash
--input <path_to_samplesheet_csv>
```

| Column   | Description |
| -------- | ----------- |
| `sample` |  Custom sample name. This identifier should be identical when you have multiple replicates from the same experimental group; just increment the `replicate` identifier appropriately. The first replicate value for any given experimental group must be `1`. Avoid including the experiment type in this identifier, since it will be parsed from the `exp_type` column and prepended to the sample name by default. |
| `fastq_1` | Full path to FastQ file for reads 1. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`. |
| `fastq_2` | Full path to FastQ file for reads 2. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`. Leave empty for single-end data. |
| `fastq_umi` | The path to the corresponding UMI `.fastq` file for deduplication. Leave empty if a separate UMI file is not available. |
| `okseq_part_file` | The path to the corresponding OK-seq partition file. Leave empty if OK-seq data is not available. Only for SCAR-seq data. |
| `replicate` | Integer representing replicate number. This will be identical for re-sequenced libraries (technical replicates). Must start from `1..<number of replicates>`. |
| `exp_type` | One of `chipseq`, `atacseq`, `scarseq`, `chorseq`. |
| `strandedness` | Either `forward` or `reverse` (SCAR-seq) or leave empty for unstranded (ChIP-seq, ChOR-seq, ATAC-seq). |
| `antibody` | This column is required to separate the downstream consensus peak merging for different antibodies. It is not advisable to generate a consensus peak set across different antibodies especially if their binding patterns are inherently different e.g. narrow transcription factors and broad histone marks. Required when `control` is specified. |
| `control` | This column should be the sample identifier for the controls for any given IP. This column together with the `control_replicate` column will set the corresponding control for each of the samples in the table. Required when `antibody` is specified. |
| `control_replicate` | Integer representing replicate number for control sample. |


### Example 1: Multiple biological replicates

This is an example of a samplesheet for a ChIP-seq experiment with one condition and two biological replicates for each antibody:

| sample | fastq_1 | fastq_2 | fastq_umi | okseq_part_file | replicate | exp_type | strandedness | antibody | control | control_replicate |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| condition_1_H3K9me3 | condition_1_bRep1_H3K9me3_R1.fastq.gz | condition_1_bRep1_H3K9me3_R3.fastq.gz | condition_1_bRep1_H3K9me3_R2.fastq.gz | | 1 | chipseq | | H3K9me3 | condition_1_INPUT | 1 |
| condition_1_H3K9me3 | condition_1_bRep2_H3K9me3_R1.fastq.gz | condition_1_bRep2_H3K9me3_R3.fastq.gz | condition_1_bRep2_H3K9me3_R2.fastq.gz | | 2 | chipseq | | H3K9me3 | condition_1_INPUT | 2 |
| condition_1_H3K27ac | condition_1_bRep1_H3K27ac_R1.fastq.gz | condition_1_bRep1_H3K27ac_R3.fastq.gz | condition_1_bRep1_H3K27ac_R2.fastq.gz | | 1 | chipseq | | H3K27ac | condition_1_INPUT | 1 |
| condition_1_H3K27ac | condition_1_bRep2_H3K27ac_R1.fastq.gz | condition_1_bRep2_H3K27ac_R3.fastq.gz | condition_1_bRep2_H3K27ac_R2.fastq.gz | | 2 | chipseq | | H3K27ac | condition_1_INPUT | 2 |
| condition_1_INPUT | condition_1_bRep1_INPUT_R1.fastq.gz | condition_1_bRep1_INPUT_R3.fastq.gz | condition_1_bRep1_INPUT_R2.fastq.gz | | 1 | chipseq | | | | |
condition_1_INPUT | condition_1_bRep2_INPUT_R1.fastq.gz | condition_1_bRep2_INPUT_R3.fastq.gz | condition_1_bRep2_INPUT_R2.fastq.gz | | 2 | chipseq | | | | |

> [!NOTE]
> You can download this example samplesheet [here](../assets/samplesheets/ex1_multiBioRep_samplesheet.csv) or copy and save the cell below:

```csv
sample,fastq_1,fastq_2,fastq_umi,okseq_part_file,replicate,exp_type,strandedness,antibody,control,control_replicate
condition_1_H3K9me3,condition_1_bRep1_H3K9me3_R1.fastq.gz,condition_1_bRep1_H3K9me3_R3.fastq.gz,condition_1_bRep1_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,condition_1_INPUT,1
condition_1_H3K9me3,condition_1_bRep2_H3K9me3_R1.fastq.gz,condition_1_bRep2_H3K9me3_R3.fastq.gz,condition_1_bRep2_H3K9me3_R2.fastq.gz,,2,chipseq,,H3K9me3,condition_1_INPUT,2
condition_1_H3K27ac,condition_1_bRep1_H3K27ac_R1.fastq.gz,condition_1_bRep1_H3K27ac_R3.fastq.gz,condition_1_bRep1_H3K27ac_R2.fastq.gz,,1,chipseq,,H3K27ac,condition_1_INPUT,1
condition_1_H3K27ac,condition_1_bRep2_H3K27ac_R1.fastq.gz,condition_1_bRep2_H3K27ac_R3.fastq.gz,condition_1_bRep2_H3K27ac_R2.fastq.gz,,2,chipseq,,H3K27ac,condition_1_INPUT,2
condition_1_INPUT,condition_1_bRep1_INPUT_R1.fastq.gz,condition_1_bRep1_INPUT_R3.fastq.gz,condition_1_bRep1_INPUT_R2.fastq.gz,,1,chipseq,,,,
condition_1_INPUT,condition_1_bRep2_INPUT_R1.fastq.gz,condition_1_bRep2_INPUT_R3.fastq.gz,condition_1_bRep2_INPUT_R2.fastq.gz,,2,chipseq,,,,
```

### Example 2: Multiple runs of the same library (technical replicates)

Both the `sample` and `replicate` identifiers have to be the same when you have sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will perform the alignments in parallel, and subsequently merge them before further analysis.

This is an example of a samplesheet for a ChIP-seq experiment with one condition, two biological replicates for each antibody, and two technical replicates for each biological replicate:

| sample | fastq_1 | fastq_2 | fastq_umi | okseq_part_file | replicate | exp_type | strandedness | antibody | control | control_replicate |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| condition_1_H3K9me3 | condition_1_bRep1_tRep1_H3K9me3_R1.fastq.gz | condition_1_bRep1_tRep1_H3K9me3_R3.fastq.gz | condition_1_bRep1_tRep1_H3K9me3_R2.fastq.gz | | 1 | chipseq | | H3K9me3 | condition_1_INPUT | 1 |
| condition_1_H3K9me3 | condition_1_bRep1_tRep2_H3K9me3_R1.fastq.gz | condition_1_bRep1_tRep2_H3K9me3_R3.fastq.gz | condition_1_bRep1_tRep2_H3K9me3_R2.fastq.gz | | 1 | chipseq | | H3K9me3 | condition_1_INPUT | 1 |
| condition_1_H3K9me3 | condition_1_bRep2_tRep1_H3K9me3_R1.fastq.gz | condition_1_bRep2_tRep1_H3K9me3_R3.fastq.gz | condition_1_bRep2_tRep1_H3K9me3_R2.fastq.gz | | 2 | chipseq | | H3K9me3 | condition_1_INPUT | 2 |
| condition_1_H3K9me3 | condition_1_bRep2_tRep2_H3K9me3_R1.fastq.gz | condition_1_bRep2_tRep2_H3K9me3_R3.fastq.gz | condition_1_bRep2_tRep2_H3K9me3_R2.fastq.gz | | 2 | chipseq | | H3K9me3 | condition_1_INPUT | 2 |
| condition_1_H3K27ac | condition_1_bRep1_tRep1_H3K27ac_R1.fastq.gz | condition_1_bRep1_tRep1_H3K27ac_R3.fastq.gz | condition_1_bRep1_tRep1_H3K27ac_R2.fastq.gz | | 1 | chipseq | | H3K27ac | condition_1_INPUT | 1 |
| condition_1_H3K27ac | condition_1_bRep1_tRep2_H3K27ac_R1.fastq.gz | condition_1_bRep1_tRep2_H3K27ac_R3.fastq.gz | condition_1_bRep1_tRep2_H3K27ac_R2.fastq.gz | | 1 | chipseq | | H3K27ac | condition_1_INPUT | 1 |
| condition_1_H3K27ac | condition_1_bRep2_tRep1_H3K27ac_R1.fastq.gz | condition_1_bRep2_tRep1_H3K27ac_R3.fastq.gz | condition_1_bRep2_tRep1_H3K27ac_R2.fastq.gz | | 2 | chipseq | | H3K27ac | condition_1_INPUT | 2 |
| condition_1_H3K27ac | condition_1_bRep2_tRep2_H3K27ac_R1.fastq.gz | condition_1_bRep2_tRep2_H3K27ac_R3.fastq.gz | condition_1_bRep2_tRep2_H3K27ac_R2.fastq.gz | | 2 | chipseq | | H3K27ac | condition_1_INPUT | 2 |
| condition_1_INPUT | condition_1_bRep1_tRep1_INPUT_R1.fastq.gz | condition_1_bRep1_tRep1_INPUT_R3.fastq.gz | condition_1_bRep1_tRep1_INPUT_R2.fastq.gz | | 1 | chipseq | | | | |
| condition_1_INPUT | condition_1_bRep1_tRep2_INPUT_R1.fastq.gz | condition_1_bRep1_tRep2_INPUT_R3.fastq.gz | condition_1_bRep1_tRep2_INPUT_R2.fastq.gz | | 1 | chipseq | | | | |
| condition_1_INPUT | condition_1_bRep2_tRep1_INPUT_R1.fastq.gz | condition_1_bRep2_tRep1_INPUT_R3.fastq.gz | condition_1_bRep2_tRep1_INPUT_R2.fastq.gz | | 2 | chipseq | | | | |
| condition_1_INPUT | condition_1_bRep2_tRep2_INPUT_R1.fastq.gz | condition_1_bRep2_tRep2_INPUT_R3.fastq.gz | condition_1_bRep2_tRep2_INPUT_R2.fastq.gz | | 2 | chipseq | | | | |

> [!NOTE]
> You can download this example samplesheet [here](../assets/samplesheets/ex2_multiTechRep_samplesheet.csv) or copy and save the cell below:

```csv
sample,fastq_1,fastq_2,fastq_umi,okseq_part_file,replicate,exp_type,strandedness,antibody,control,control_replicate
condition_1_H3K9me3,condition_1_bRep1_tRep1_H3K9me3_R1.fastq.gz,condition_1_bRep1_tRep1_H3K9me3_R3.fastq.gz,condition_1_bRep1_tRep1_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,condition_1_INPUT,1
condition_1_H3K9me3,condition_1_bRep1_tRep2_H3K9me3_R1.fastq.gz,condition_1_bRep1_tRep2_H3K9me3_R3.fastq.gz,condition_1_bRep1_tRep2_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,condition_1_INPUT,1
condition_1_H3K9me3,condition_1_bRep2_tRep1_H3K9me3_R1.fastq.gz,condition_1_bRep2_tRep1_H3K9me3_R3.fastq.gz,condition_1_bRep2_tRep1_H3K9me3_R2.fastq.gz,,2,chipseq,,H3K9me3,condition_1_INPUT,2
condition_1_H3K9me3,condition_1_bRep2_tRep2_H3K9me3_R1.fastq.gz,condition_1_bRep2_tRep2_H3K9me3_R3.fastq.gz,condition_1_bRep2_tRep2_H3K9me3_R2.fastq.gz,,2,chipseq,,H3K9me3,condition_1_INPUT,2
condition_1_H3K27ac,condition_1_bRep1_tRep1_H3K27ac_R1.fastq.gz,condition_1_bRep1_tRep1_H3K27ac_R3.fastq.gz,condition_1_bRep1_tRep1_H3K27ac_R2.fastq.gz,,1,chipseq,,H3K27ac,condition_1_INPUT,1
condition_1_H3K27ac,condition_1_bRep1_tRep2_H3K27ac_R1.fastq.gz,condition_1_bRep1_tRep2_H3K27ac_R3.fastq.gz,condition_1_bRep1_tRep2_H3K27ac_R2.fastq.gz,,1,chipseq,,H3K27ac,condition_1_INPUT,1
condition_1_H3K27ac,condition_1_bRep2_tRep1_H3K27ac_R1.fastq.gz,condition_1_bRep2_tRep1_H3K27ac_R3.fastq.gz,condition_1_bRep2_tRep1_H3K27ac_R2.fastq.gz,,2,chipseq,,H3K27ac,condition_1_INPUT,2
condition_1_H3K27ac,condition_1_bRep2_tRep2_H3K27ac_R1.fastq.gz,condition_1_bRep2_tRep2_H3K27ac_R3.fastq.gz,condition_1_bRep2_tRep2_H3K27ac_R2.fastq.gz,,2,chipseq,,H3K27ac,condition_1_INPUT,2
condition_1_INPUT,condition_1_bRep1_tRep1_INPUT_R1.fastq.gz,condition_1_bRep1_tRep1_INPUT_R3.fastq.gz,condition_1_bRep1_tRep1_INPUT_R2.fastq.gz,,1,chipseq,,,,
condition_1_INPUT,condition_1_bRep1_tRep2_INPUT_R1.fastq.gz,condition_1_bRep1_tRep2_INPUT_R3.fastq.gz,condition_1_bRep1_tRep2_INPUT_R2.fastq.gz,,1,chipseq,,,,
condition_1_INPUT,condition_1_bRep2_tRep1_INPUT_R1.fastq.gz,condition_1_bRep2_tRep1_INPUT_R3.fastq.gz,condition_1_bRep2_tRep1_INPUT_R2.fastq.gz,,2,chipseq,,,,
condition_1_INPUT,condition_1_bRep2_tRep2_INPUT_R1.fastq.gz,condition_1_bRep2_tRep2_INPUT_R3.fastq.gz,condition_1_bRep2_tRep2_INPUT_R2.fastq.gz,,2,chipseq,,,,
```

### Example 3: Full design

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet.

A final design file may look something like the one below. This is for two antibodies and associated controls, where the second replicate of the `WT_BCATENIN_IP` and `NAIVE_BCATENIN_IP` samples have been sequenced twice:

> [!NOTE]
> You can download this example samplesheet [here]() or copy and save the cell below:

## Reference genome files

### Genome index, FASTA, GTF/GFF and gene BED files

The minimum reference genome requirements are a FASTA and a GTF file, all other files required to run the pipeline can be generated from these files. However, it is more storage and compute friendly if you are able to re-use reference genome files as efficiently as possible. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices (e.g. those unavailable on [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space. You can then either provide the appropriate reference genome files on the command-line via the appropriate parameters (e.g. `--bwa_index '/path/to/bwa/index/'`) or via a [custom config file](https://nf-co.re/usage/configuration#custom-configuration-files).

<!-- - If `--genome` is provided then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`. -->
- If `--gene_bed` is not provided then it will be generated from the GTF file.

> **NB:** Compressed reference files are also supported by the pipeline i.e. standard files with the `.gz` extension and indices folders with the `tar.gz` extension.

### Genome blacklist regions

The blacklist bed files where obtained using the commands below:

```console
cd ..
mkdir -p v1.0
cd v1.0
wget -L https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz && gunzip ENCFF001TDO.bed.gz && mv ENCFF001TDO.bed hg19-blacklist.v1.bed

mkdir -p assets/blacklists/v2.0/
cd assets/blacklists/v2.0/
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.v2.bed.gz && gunzip ce10-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.v2.bed.gz && gunzip ce11-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.v2.bed.gz && gunzip dm3-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz && gunzip dm6-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz && gunzip hg19-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz && gunzip hg38-blacklist.v2.bed.gz
wget -L https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz && gunzip mm10-blacklist.v2.bed.gz

cd ..
mkdir -p v3.0
cd v3.0
wget -L https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz && gunzip ENCFF356LFX.bed.gz && mv ENCFF356LFX.bed hg38-blacklist.v3.bed
```

> **NB:** A detailed description of the different versions of the files can be found [here](https://sites.google.com/site/anshulkundaje/projects/blacklists). Also, to to see which blacklist bed files are assigned by default to the respective reference genome check the [igenomes.config](conf/igenomes.config).


### Initiation zones
-

### OK-seq partitions
-

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run grothlab/glseq \
  -profile docker \
  --input samplesheet.csv \
  --genome GRCh37 \
  --outdir <OUTDIR>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
> The above pipeline run specified with a params file in yaml format:

```bash
nextflow run grothlab/glseq -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
input: 'data'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Parameters

#### Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 11 columns, and a header row. See [usage docs](usage.md).</small></details>| `string` |  |  |  |
| `fragment_size` | Estimated fragment size used to extend single-end reads. | `integer` | 150 |  |  |
| `seq_platform` | Platform/technology used to produce the reads. Corresponds to the `PL` tag in the SAM/BAM file header. <details><summary>Help</summary><small>See the [SAM format specification](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf). Valid values: CAPILLARY, DNBSEQ (MGI/BGI), ELEMENT, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT (Oxford Nanopore), PACBIO (Pacific Biosciences), SINGULAR, SOLID, and ULTIMA. This field should be omitted when the technology is not in this list (though the PM field may still be present in this case) or is unknown.</small></details>| `string` | None |  |  |
| `seq_center` | Sequencing center information to be added to read group of BAM files. | `string` |  |  |  |
| `read_length` | Read length used to calculate MACS3 genome size for peak calling if `--macs_gsize` isn't provided. | `integer` | 50 |  |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

#### Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `genome` | Name of iGenomes reference. <details><summary>Help</summary><small>If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`. <br><br>See the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details.</small></details>| `string` |  |  |  |
| `fasta` | Path to FASTA genome file. <details><summary>Help</summary><small>This parameter is *mandatory* if `--genome` is not specified. If you don't have the appropriate alignment index available this will be generated for you automatically. Combine with `--save_reference` to save alignment index for future runs.</small></details>| `string` |  |  |  |
| `gtf` | Path to GTF annotation file. <details><summary>Help</summary><small>This parameter is *mandatory* if `--genome` is not specified.</small></details>| `string` |  |  |  |
| `gff` | Path to GFF3 annotation file. <details><summary>Help</summary><small>This parameter must be specified if `--genome` or `--gtf` are not specified.</small></details>| `string` |  |  |  |
| `bwa_index` | Path to directory or tar.gz archive for pre-built BWA index. | `string` |  |  |  |
| `bowtie2_index` | Path to directory or tar.gz archive for pre-built Bowtie2 index. | `string` |  |  |  |
| `chromap_index` | Path to directory or tar.gz archive for pre-built Chromap index. | `string` |  |  |  |
| `star_index` | Path to directory or tar.gz archive for pre-built STAR index. | `string` |  |  |  |
| `gene_bed` | Path to BED file containing gene intervals. This will be created from the GTF file if not specified. | `string` |  |  |  |
| `macs_gsize` | Effective genome size parameter required by MACS3. <details><summary>Help</summary><small>[Effective genome size](https://github.com/taoliu/MACS#-g--gsize) parameter required by MACS3. If using an iGenomes reference these have been provided when `--genome` is set as *GRCh37*, *GRCh38*, *GRCm38*, *WBcel235*, *BDGP6*, *R64-1-1*, *EF2*, *hg38*, *hg19* and *mm10*. For other genomes, if this parameter is not specified then the MACS3 peak-calling and differential analysis will be skipped.</small></details>| `number` |  |  |  |
| `blacklist` | Path to blacklist regions in BED format, used for filtering alignments. <details><summary>Help</summary><small>If provided, alignments that overlap with the regions in this file will be filtered out (see [ENCODE blacklists](https://sites.google.com/site/anshulkundaje/projects/blacklists)). The file should be in BED format. Blacklisted regions for *GRCh37*, *GRCh38*, *GRCm38*, *hg19*, *hg38*, *mm10* are bundled with the pipeline in the [`blacklists`](../assets/blacklists/) directory, and as such will be automatically used if any of those genomes are specified with the `--genome` parameter.</small></details>| `string` |  |  |  |
| `save_reference` | If generated by the pipeline save the BWA index in the results directory. <details><summary>Help</summary><small>If the BWA index is generated by the pipeline use this parameter to save it to your results folder. These can then be used for future pipeline runs, reducing processing times.</small></details>| `boolean` |  |  |  |
| `igenomes_base` | Directory / URL base for iGenomes references. | `string` | s3://ngi-igenomes/igenomes/ |  | True |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>| `boolean` | True |  | True |

#### UMI extraction and deduplication options

Options to adjust UMI extraction and deduplication criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_umi` | Enable UMI extraction and deduplication. <details><summary>Help</summary><small>Set this parameter to enable UMI extraction and deduplication. </small></details>| `boolean` |  |  |  |
| `skip_umi_extract` | Skip UMI extraction step. <details><summary>Help</summary><small>Use this flag to skip the UMI extraction/transfer step. This is useful if you have already extracted UMIs from your data but want to run the deduplication</small></details>| `boolean` |  |  |  |
| `umi_discard_read` | Discard R1 or R2 if required | `integer` |  |  |  |
| `get_dedup_stats` | Get deduplication statistics. | `boolean` | True |  |  |

#### Adapter trimming options

Options to adjust adapter trimming criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `clip_r1` | Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads). | `integer` |  |  |  |
| `clip_r2` | Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only). | `integer` |  |  |  |
| `three_prime_clip_r1` | Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed. | `integer` |  |  |  |
| `three_prime_clip_r2` | Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed. | `integer` |  |  |  |
| `trim_nextseq` | Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails. <details><summary>Help</summary><small>This enables the option Cutadapt `--nextseq-trim=3'CUTOFF` option via Trim Galore, which will set a quality cutoff (that is normally given with -q instead), but qualities of G bases are ignored. This trimming is in common for the NextSeq- and NovaSeq-platforms, where basecalls without any signal are called as high-quality G bases.</small></details>| `integer` |  |  |  |
| `skip_trimming` | Skip the adapter trimming step. <details><summary>Help</summary><small>Use this if your input FastQ files have already been trimmed outside of the workflow or if you're very confident that there is no adapter contamination in your data.</small></details>| `boolean` |  |  |  |
| `save_trimmed` | Save the trimmed FastQ files in the results directory. <details><summary>Help</summary><small>By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.</small></details>| `boolean` |  |  |  |
| `min_trimmed_reads` | Minimum number of reads required after trimming. | `integer` | 1 |  |  |

#### Alignment options

Options to adjust parameters and filtering criteria for read alignments.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `aligner` | Specifies the alignment algorithm to use - available options are 'bwa', 'bowtie2' and 'star'. | `string` | bowtie2 |  |  |
| `bwa_min_score` | Don’t output BWA MEM alignments with score lower than this parameter. | `integer` |  |  |  |
| `sort_bam` | Should the BAM files be sorted by coordinate? Only relevant if Bowtie2 or BWA are used as aligners. | `boolean` | True |  |  |
| `save_align_intermeds` | Save the intermediate BAM files from the alignment step. <details><summary>Help</summary><small>By default, intermediate BAM files will not be saved. The final BAM files created after the appropriate filtering step are always saved to limit storage usage. Set this parameter to also save other intermediate BAM files.</small></details>| `boolean` |  |  |  |
| `save_unaligned` | Where possible, save unaligned reads from either STAR, HISAT2 or Salmon to the results directory. <details><summary>Help</summary><small>This may either be in the form of FastQ or BAM files depending on the options available for that particular tool.</small></details>| `boolean` |  |  |  |

#### Spike-in normalization options

Options to adjust spike-in splitting and normalization criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `spikein_genome` | Name of the spike-in genome if used in the experiment. The spike-in chromosome names in the `--fasta` and `<aligner>_index` files should have the format `<chromosome_name>_<spikein_genome>`, e.g. `chr1_dm6`. | `string` | None |  |  |
| `save_spikein_intermeds` | Save intermediate files for spike-in genome alignment. | `boolean` | False |  |  |

#### Multimapper allocation options

Options to adjust multimapper allocation criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `allocate_n_multimappers` | Max number of multimappers to output from alignment to allocate. | `integer` | 0 |  |  |
| `allocation_method` | Method to allocate multimappers. | `string` | allo |  |  |
| `allo_mixed_cnn` | Use CNN trained on histone ChIP-seq datasets with mixed peaks, narrow by default. See [allo](https://github.com/seqcode/allo?tab=readme-ov-file#options) for more details. | `boolean` | True |  |  |
| `save_allocation_intermeds` | Save the intermediate BAM files from the multimapper allocation step. | `boolean` | True |  |  |

#### Peak calling options

Options to adjust peak calling criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `narrow_peak` | Run MACS3 in narrowPeak mode. <details><summary>Help</summary><small>MACS3 is run by default with the [`--broad`](https://github.com/taoliu/MACS#--broad) flag. Specify this flag to call peaks in narrowPeak mode.</small></details>| `boolean` |  |  |  |
| `broad_cutoff` | Specifies broad cutoff value for MACS3. Only used when --narrow_peak isnt specified. | `number` | 0.1 |  |  |
| `macs_fdr` | Minimum FDR (q-value) cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive. | `number` |  |  |  |
| `macs_pvalue` | p-value cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive. If --macs_pvalue cutoff is set, q-value will not be calculated and reported as -1 in the final .xls file. | `number` |  |  |  |
| `min_reps_consensus` | Number of biological replicates required from a given condition for a peak to contribute to a consensus peak. <details><summary>Help</summary><small>If you are confident you have good reproducibility amongst your replicates then you can increase the value of this parameter to create a 'reproducible' set of consensus peaks. For example, a value of 2 will mean peaks that have been called in at least 2 replicates will contribute to the consensus set of peaks, and as such peaks that are unique to a given replicate will be discarded.</small></details>| `integer` | 1 |  |  |
| `save_macs_pileup` | Instruct MACS3 to create bedGraph files normalised to signal per million reads. | `boolean` | True |  |  |
| `skip_peak_qc` | Skip MACS3 peak QC plot generation. | `boolean` |  |  |  |
| `skip_peak_annotation` | Skip annotation of MACS3 and consensus peaks with HOMER. | `boolean` |  |  |  |
| `skip_consensus_peaks` | Skip consensus peak generation, annotation and counting. | `boolean` |  |  |  |
| `save_macs_pileup_intermeds` |  | `boolean` |  |  |  |

#### Process skipping options

Options to skip various steps within the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_fastqc` | Skip FastQC. | `boolean` |  |  |  |
| `skip_picard_metrics` | Skip Picard CollectMultipleMetrics. | `boolean` |  |  |  |
| `skip_preseq` | Skip Preseq. | `boolean` | True |  |  |
| `deseq2_vst` | Use vst transformation instead of rlog with DESeq2. <details><summary>Help</summary><small>See [DESeq2 docs](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization).</small></details>| `boolean` | True |  |  |
| `skip_plot_profile` | Skip deepTools plotProfile. | `boolean` |  |  |  |
| `skip_plot_fingerprint` | Skip deepTools plotFingerprint. | `boolean` |  |  |  |
| `skip_spp` | Skip Phantompeakqualtools. | `boolean` |  |  |  |
| `skip_deseq2_qc` | Skip DESeq2 PCA and heatmap plotting. | `boolean` |  |  |  |
| `skip_igv` | Skip IGV. | `boolean` |  |  |  |
| `skip_multiqc` | Skip MultiQC. | `boolean` |  |  |  |
| `skip_qc` | Skip all QC steps except for MultiQC. | `boolean` |  |  |  |

#### SCAR-seq analysis options

Options to adjust SCAR-seq analysis criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `scar_slop` |  | `integer` | 0 |  |  |
| `scar_window_size` |  | `integer` | 1000 |  |  |
| `scar_step_size` |  | `integer` | 1000 |  |  |
| `scar_radius` |  | `integer` | 30 |  |  |
| `scar_dradius` |  | `integer` | 30 |  |  |
| `scar_zradius` |  | `integer` | 1 |  |  |
| `scar_cpm_cutoff` |  | `number` | 0.3 |  |  |
| `scar_plot_range` |  | `integer` | 100 |  |  |
| `scar_exclude_chromosomes` |  | `string` | chrX,chrY,chrM |  |  |
| `save_scarseq_intermeds` |  | `boolean` |  |  |  |

#### Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

#### Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  | True |

#### Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `fingerprint_bins` | Number of genomic bins to use when calculating deepTools fingerprint plot. | `integer` | 500000 |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |

#### Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `multiqc_logo` |  | `string` |  |  |  |
| `multiqc_methods_description` |  | `string` |  |  |  |
| `hook_url` |  | `string` |  |  |  |
| `version` |  | `boolean` |  |  |  |
| `pipelines_testdata_base_path` |  | `string` | https://raw.githubusercontent.com/nf-core/test-datasets/ |  |  |
| `local_testdata_base_path` |  | `string` | /maps/projects/dan1/data/Groth_group/SRP/glseq_testdata/ |  |  |
| `validationFailUnrecognisedParams` |  | `boolean` |  |  |  |
| `validationLenientMode` |  | `boolean` |  |  |  |
| `validationShowHiddenParams` |  | `boolean` |  |  |  |


### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull grothlab/glseq
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

<!-- First, go to the [nf-core/chipseq releases page](https://github.com/nf-core/chipseq/releases) and find the latest pipeline version - numeric only (eg. `2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0.0`. Of course, you can switch to another version by changing the number after the `-r` flag. -->

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

### Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

### Custom configuration

#### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

#### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

#### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

#### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

#### Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```


## Quick start guide for DAN System users

1. Read the [DAN System User Guide](https://sgn102.pages.ku.dk/a-not-long-tour-of-dangpu/) to understand how to use the DAN System.

2. Start a [*tmux*](https://github.com/tmux/tmux/wiki/Getting-Started) session:

    ```bash
    tmux new-session -s <session_name>
    ```

3. Launch a minimal interactive [*slurm*](https://slurm.schedmd.com/documentation.html) job session:

    ```bash
    srun -c 1 --mem=2gb --time=2-00:00:00 --pty bash
    ```
> [!NOTE]  
> Adjust `--time` as necessary, the command above keeps the slurm job active for two days (more than enough for most pipeline runs).

4. Load the required [*modules*](https://modules.readthedocs.io/en/latest/):

    ```bash
    module load openjdk/20.0.0 nextflow/24.04.4 singularity/3.8.7
    ```

5. Clone the pipeline repository:

  - Generate a Personal Access Token (PAT):

      - Go to [GitHub](https://github.com/) and log in to your account.

      - Make sure you are part of the [Groth Lab organization](https://github.com/grothlab). If not, please contact Nicolás Alcaraz ([nicolas.alcaraz@cpr.ku.dk](nicolas.alcaraz@cpr.ku.dk)) to request access.

      - Click on the profile picture in the right-hand menu, then ***Settings*** > ***Developer settings*** > [***Tokens (classic)***](https://github.com/settings/tokens).

      - Click on the ***Generate new token*** and ***Generate new token (classic)*** buttons.

      - Provide a meaningful name to identify its purpose (e.g. `glseq_pat`) and select the required permissions: for cloning glseq, the “repo” permissions are sufficient.

      - Click on the ***Generate Token*** button to generate your PAT.
      - Copy the generated token to your clipboard. Remember that PATs are sensitive and should be treated like passwords.

>[!WARNING]
> Make note of the token because once you close the window you won’t be able to view the token again!

  - Open your terminal and create or navigate to the directory where you want to clone the repository:

      ```bash
      mkdir -p <path_to_software_directory>
      cd <path_to_software_directory>
      ```

  - Clone the repository using the personal access token:

      ```bash
      git clone https://github.com/grothlab/glseq.git
      ```
    
      ```bash
      username : <your_username>
      password : <your_generated_token>
      ```
    
6. Create an output directory for your pipeline run if it does not exist:

    ```bash
    mkdir -p <path_to_output_directory>
    cd <path_to_output_directory>
    ```

7. Run a pipeline test (`local_test_scarseq`, `local_test_chipseq`, `local_test_atacseq` or `local_test_chorseq`) with the institution profile ([`ku_sund_danhead_mod`](conf/ku_sund_danhead_mod.config)):

    ```bash
    nextflow run <path_to_software_directory>/glseq \
      -profile ku_sund_danhead_mod,local_test_scarseq \
      --outdir <path_to_output_directory>
    ```

  > [!TIP]
  >  Include the `-work-dir` argument if you want to save the work/temporary files in a specific directory to inspect them later. Otherwise, these files are saved in `/scratch/temp/$::env(USER)/nxf/work` by default.

  > [!NOTE]
  > The [`ku_sund_danhead`](https://github.com/nf-core/configs/blob/master/docs/ku_sund_danhead.md) [config profile](https://github.com/nf-core/configs/blob/master/conf/ku_sund_danhead.config) created by the DAN System administrators has set up [`cleanup = true`](https://www.nextflow.io/docs/stable/reference/config.html#unscoped-options) by default, which automatically deletes all files in the work directory on a "successful" completion of a run. However, this prevents the use of the `-resume` feature on subsequent executions of any pipeline run, and thus, a modified version ([`ku_sund_danhead_mod`](conf/ku_sund_danhead_mod.config)) with `cleanup = false` was created to facilitate the running and resuming of this pipeline.
  

8. You can now detach from the *tmux* session by pressing `Ctrl+b` and then `d`. You can reattach to the session later by running:

    ```bash
    tmux attach-session -t <session_name>
    ```

9. Run your own analysis, for example:

    ```bash
    nextflow run <path_to_software_directory>/glseq \
      -profile ku_sund_danhead_mod \
      --input <path_to_your_input_samplesheet_csv_file> \
      --with_umi \
      --skip_umi_extract false \
      --genome mm10 \
      --spikein_genome dm6 \
      ...
      --outdir <path_to_output_directory>
    ```

### Reference genome files

Reference files for the Groth Lab have been made available by [Nicolás Alcaraz](https://github.com/satroz). If you have any questions, requests, or need help with the reference files, please check these [slides](https://alumni-my.sharepoint.com/:p:/g/personal/ngl887_ku_dk/EXx9wPQ1shlNtnVoafTB5iQBnyYvGUZZf3ckArr0cZqiIQ?e=xQEgKb) or contact him ([nicolas.alcaraz@cpr.ku.dk](nicolas.alcaraz@cpr.ku.dk)).

#### Genome index, FASTA, GTF/GFF and gene BED files

- Prebuilt indexes are inputted using the `--<aligner>_index` parameter, and can be found at:

  ```
  /maps/projects/dan1/data/Groth_group/shared/references/<organism>/<genome_version>/indices/<aligner>/
  ```

- Genome FASTA files are inputted using the `--fasta` parameter, and can be found at:

  ```
  /maps/projects/dan1/data/Groth_group/shared/references/<organism>/<genome_version>/genome/fasta/*.fa.gz
  ```

- GTF/GFF files are inputted using the `--gtf` or `--gff` parameters, and can be found at:

  ```
  /maps/projects/dan1/data/Groth_group/shared/references/<organism>/<genome_version>/annotations/transcript_models/<gtf_or_gff>/*.<gtf_or_gff>.gz
  ```

- Gene BED files are inputted using the `--gene_bed` parameter, and can be found at:

  ```
  /maps/projects/dan1/data/Groth_group/shared/references/<organism>/<genome_version>/annotations/transcript_models/*.bed.gz
  ```

> [!TIP]
> Make sure to always input the prebuilt index files (besides the FASTA and GTF files) if available to avoid building the index from scratch every time you run the pipeline:
> ```bash
> nextflow run grothlab/glseq \
>   --bowtie2_index <path_to_index> \
>   --fasta <path_to_fasta> \
>   --gtf <path_to_gtf>
>```


#### Genome blacklist regions

Blacklist BED files are inputted using the `--blacklist` parameter, and can be found at:

```
/maps/projects/dan1/data/Groth_group/shared/references/<organism>/<genome_version>/annotations/blacklists/bed/*.bed.gz
```

#### Initiation zones

Initiation zone BED files are inputted using the `--initiation_zones` parameter, and can be found at:

```bash
/maps/projects/dan1/data/Groth_group/shared/references/<organism>/<genome_version>/external_data/Replication/Okasaki_seq/bed_files/OKseq_Initiation_Zones_*.bed.gz
```

#### OK-seq partitions

Okazaki fragment sequencing (OK-seq) partition BED files are inputted for each SCAR-seq sample through the [`okseq_part_file` column](#samplesheet-input) in the samplesheet, and can be found at:

```bash
/maps/projects/dan1/data/Groth_group/shared/references/<organism>/<genome_version>/external_data/Replication/Okasaki_seq/rfd_files/OKseq_RFD_*.bed.gz
```

#### Example run using local genome files

On the command line:

```bash
nextflow run /user/datadir/software/glseq \
      -profile ku_sund_danhead_mod \
      --input /user/datadir/projects/project1/project1_glseq_samplesheet.csv \
      --outdir /user/datadir/projects/project1/output/ \
      --work-dir /user/datadir/projects/project1/output/work/ \
      --aligner bowtie2 \
      --read_length 50 \
      --with_umi \
      --skip_umi_extract false \
      --genome mm10 \
      --spikein_genome dm6 \
      --bowtie2_index /home/rlh546/Groth_group/shared/references/Mus_musculus/GRCm38/spiked/GRCm38_dm6/indices/bowtie2_2.5.4_index \
      --fasta /maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/spiked/GRCm38_dm6/genome/fasta/GRCm38_dm6.primary_assembly.genome.fa.gz \
      --gtf /maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/annotations/transcript_models/gtf/gencode.vM25.primary_assembly.annotation.gtf.gz \
      --initiation_zones /maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/external_data/Replication/Okazaki_seq/bed_files/OKseq_Initiation_Zones_mESC_SRR7535256_r1_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz \
      --save_align_intermeds \
      --save_spikein_intermeds
```

### Tips

Parameters specified on the command line [can be also specified in a params file using the `-params-file` option](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters). This is useful for saving a set of parameters that you use frequently, or for sharing parameters with others.


Parameters can be represented in YAML (`.yml`) format:

```bash
nextflow run /user/datadir/software/glseq \
      -profile ku_sund_danhead_mod \
      --params-file /user/datadir/projects/project1/project1_glseq_params.yml
```

```yaml title="project1_glseq_params.yml"
input: /user/datadir/projects/project1/project1_glseq_samplesheet.csv
outdir: /user/datadir/projects/project1/output/
work_dir: /user/datadir/projects/project1/output/work/
aligner: bowtie2
read_length: 50
with_umi: true
skip_umi_extract: false
genome: mm10
spikein_genome: dm6
bowtie2_index: /home/rlh546/Groth_group/shared/references/Mus_musculus/GRCm38/spiked/GRCm38_dm6/indices/bowtie2_2.5.4_index
fasta: /maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/spiked/GRCm38_dm6/genome/fasta/GRCm38_dm6.primary_assembly.genome.fa.gz
gtf: /maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/annotations/transcript_models/gtf/gencode.vM25.primary_assembly.annotation.gtf.gz
initiation_zones: /maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/external_data/Replication/Okazaki_seq/bed_files/OKseq_Initiation_Zones_mESC_SRR7535256_r1_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz
save_align_intermeds: true
save_spikein_intermeds: true
```

Or in JSON (`.json`) format:

```bash
nextflow run /user/datadir/software/glseq \
      -profile ku_sund_danhead_mod \
      --params-file /user/datadir/projects/project1/project1_glseq_params.json
```

```json title="project1_glseq_params.json"
{
  "input": "/user/datadir/projects/project1/project1_glseq_samplesheet.csv",
  "outdir": "/user/datadir/projects/project1/output/",
  "work_dir": "/user/datadir/projects/project1/output/work/",
  "aligner": "bowtie2",
  "read_length": 50,
  "with_umi": true,
  "skip_umi_extract": false,
  "genome": "mm10",
  "spikein_genome": "dm6",
  "bowtie2_index": "/home/rlh546/Groth_group/shared/references/Mus_musculus/GRCm38/spiked/GRCm38_dm6/indices/bowtie2_2.5.4_index",
  "fasta": "/maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/spiked/GRCm38_dm6/genome/fasta/GRCm38_dm6.primary_assembly.genome.fa.gz",
  "gtf": "/maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/annotations/transcript_models/gtf/gencode.vM25.primary_assembly.annotation.gtf.gz",
  "initiation_zones": "/maps/projects/dan1/data/Groth_group/shared/references/Mus_musculus/GRCm38/external_data/Replication/Okazaki_seq/bed_files/OKseq_Initiation_Zones_mESC_SRR7535256_r1_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz",
  "save_align_intermeds": true,
  "save_spikein_intermeds": true
}
```