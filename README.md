[![GitHub Actions CI Status](https://github.com/grothlab/glseq/actions/workflows/ci.yml/badge.svg)](https://github.com/grothlab/glseq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/grothlab/glseq/actions/workflows/linting.yml/badge.svg)](https://github.com/grothlab/glseq/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/grothlab/glseq)

## Introduction

**grothlab/glseq** is a bioinformatics pipeline for the analysis of sequencing data (ChIP-seq, [SCAR-seq](https://doi.org/10.1038/s41596-021-00585-3), [ChOR-seq](https://doi.org/10.1038/s41596-021-00585-3), etc.).

<!-- On release, automated continuous integration tests run the pipeline on a [full-sized dataset](https://github.com/nf-core/test-datasets/tree/chipseq#full-test-dataset-origin) on the AWS cloud infrastructure. The dataset consists of FoxA1 (transcription factor) and EZH2 (histone,mark) IP experiments from _Franco et al. 2015_ ([GEO: GSE59530](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59530), [PMID: 25752574](https://pubmed.ncbi.nlm.nih.gov/25752574/)) and _Popovic et al. 2014_ ([GEO: GSE57632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57632), [PMID: 25188243](https://pubmed.ncbi.nlm.nih.gov/25188243/)), respectively. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from running the full-sized tests can be viewed on the [nf-core website](https://nf-co.re/chipseq/results). -->

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. [Optional] UMI extraction ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools) or [`umi-transfer`](https://github.com/SciLifeLab/umi-transfer))
3. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
4. Choice of multiple aligners:
    1. ([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
    2. ([`Chromap`](https://github.com/haowenz/chromap))
    3. ([`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
    4. ([`STAR`](https://github.com/alexdobin/STAR))
5. Merge alignments from multiple libraries of the same sample (technical replicates) ([`picard`](https://broadinstitute.github.io/picard/))
6. [Optional] Marking duplicates ([`picard`](https://broadinstitute.github.io/picard/))


7. BAM filtering to remove:
    - reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
    - reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    - reads that are not marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    - reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    - reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    - reads containing > 4 mismatches ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    - reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools); _paired-end only_)
    - reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
    - reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
    - reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
8. Alignment-level QC and estimation of library complexity ([`picard`](https://broadinstitute.github.io/picard/), [`Preseq`](http://smithlabresearch.org/software/preseq/))
9. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
10. Generate gene-body meta-profile from bigWig files ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html))
11. Calculate genome-wide IP enrichment relative to control ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html))
12. Calculate strand cross-correlation peak and ChIP-seq quality measures including NSC and RSC ([`phantompeakqualtools`](https://github.com/kundajelab/phantompeakqualtools))

13. ChIP-seq downstream analysis:

    1. Call broad/narrow peaks ([`MACS3`](https://github.com/macs3-project/MACS))
    2. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
    3. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
    4. Count reads in consensus peaks ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
    5. PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    6. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
    7. Present QC for raw read, alignment, peak-calling and differential binding results ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

14. SCAR-seq downstream analysis:

    1. Splitting BAM files by forward and reverse strands.
    2. Computing BEDGRAPH summaries of feature coverage per strand ([`BEDTools`](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html))
    3. Creating BigWig files from BEDGRAPH files ([`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    4. Splitting chromosome windows
    5. Computing the average coverage per window ([`bigWigAverageOverBed`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    6. Normalizing per strand and chromosome with counts per million (CPM)
    7. Calculating replication fork directionality (RFD) ([`partition_smooth.pl`]())
    8. Generating partition files for samples, stranded inputs and input-corrected samples.
    9. Plotting partition files and scatter-correlation plots against OK-seq if provided.


## Quick start for DAN System users

1. Read the [DAN System User Guide](https://sgn102.pages.ku.dk/a-not-long-tour-of-dangpu/) to understand how to use the DAN System.

2. Start a [*tmux*](https://github.com/tmux/tmux/wiki/Getting-Started) session:

    ```bash
    tmux new-session -s <session_name>
    ```

3. Launch a minimal interactive [*slurm*](https://slurm.schedmd.com/documentation.html) job session:

    ```bash
    srun -c 1 --mem=1gb --time=6-00:00:00 --pty bash
    ```

4. Load the required [*modules*](https://modules.readthedocs.io/en/latest/):

    ```bash
    module load openjdk/20.0.0 nextflow/24.04.4 singularity/3.8.7
    ```
5. Create output directory if it does not exist:

    ```bash
    mkdir -p <path_to_output_directory>
    cd <path_to_output_directory>
    ```

6. Run a pipeline test (`local_test_scarseq`, `local_test_chipseq`, or `local_test_full`) with the institution profile ([`ku_sund_danhead`](https://github.com/nf-core/configs/blob/master/docs/ku_sund_danhead.md)):

    ```bash
    nextflow run grothlab/glseq \
      -profile ku_sund_danhead,local_test_scarseq \
      --outdir <path_to_output_directory>
    ```

7. You can now detach from the *tmux* session by pressing `Ctrl+b` and then `d`. You can reattach to the session later by running:

    ```bash
    tmux attach-session -t <session_name>
    ```

8. Run your own analysis:

    ```bash
    nextflow run grothlab/glseq \
      -profile ku_sund_danhead \
      --input <path_to_your_input_samplesheet.csv> \
      --outdir <path_to_output_directory>
    ```

## Quick start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run grothlab/glseq \
      -profile test,YOURPROFILE \
      --outdir <path-of-output-directory>
    ```

    Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

    > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
    > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
    > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    ```bash
    nextflow run grothlab/glseq \
      --input <path-of-your-input-samplesheet-csv-file> \
      --outdir <path-of-output-directory> \
      --genome GRCh37 \
      -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

See [usage docs](https://nf-co.re/chipseq/usage) for all of the available options when running the pipeline.

## Usage

### Samplesheet input



You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the `--input` parameter to specify its location. It has to be a comma-separated (`.csv`) file with with 11 columns and a header row as explained below.

| Column   | Description |
| -------- | ----------- |
| `sample` |  This identifier should be identical when you have multiple replicates from the same experimental group; just increment the `replicate` identifier appropriately. The first replicate value for any given experimental group must be `1`. |
| `fastq_1` | Full path to FastQ file for reads 1. File has to be gzipped and have the extension “.fastq.gz” or “.fq.gz”. |
| `fastq_2` | Full path to FastQ file for reads 2. File has to be gzipped and have the extension “.fastq.gz” or “.fq.gz”. Leave empty for single-end data. |
| `fastq_umi` | The path to the corresponding UMI `.fastq` file for deduplication. Leave empty if a separate UMI file is not available. |
| `okseq_part_file` | The path to the corresponding OK-seq partition file. Leave empty if OK-seq data is not available. Only for SCAR-seq data. |
| `replicate` | Integer representing replicate number. This will be identical for re-sequenced libraries. Must start from `1..<number of replicates>`. |
| `exp_type` | Either `chipseq` or `scarseq`. |
| `strandedness` | Either `forward`, `reverse` or leave empty for unstranded (ChIP-seq). |
| `antibody` | This column is required to separate the downstream consensus peak merging for different antibodies. It is not advisable to generate a consensus peak set across different antibodies especially if their binding patterns are inherently different e.g. narrow transcription factors and broad histone marks. Required when control is specified. |
| `control` | This column should be the sample identifier for the controls for any given IP. This column together with the `control_replicate` column will set the corresponding control for each of the samples in the table. |
| `control_replicate` | Integer representing replicate number for control sample. |

<details open>
<summary>
<b>Here is an example samplesheet for running the pipeline:</b>

**[samplesheet_template.csv](https://github.com/grothlab/glseq/blob/dev/assets/samplesheet_template_scarseq.csv):**

```csv
sample,fastq_1,fastq_2,fastq_umi,okseq_part_file,replicate,exp_type,strandedness,antibody,control,control_replicate
project1_scar_cond1_H3K9me3,/path/to/samples/project1_scar_cond1_H3K9me3_R1.fastq.gz,/path/to/samples/project1_scar_cond1_H3K9me3_R3.fastq.gz,/path/to/samples/project1_scar_cond1_H3K9me3_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,H3K9me3,project1_scar_cond1_strandedInput,1
project1_scar_cond1_H4K20me0,/path/to/samples/project1_scar_cond1_H4K20me0_R1.fastq.gz,/path/to/samples/project1_scar_cond1_H4K20me0_R3.fastq.gz,/path/to/samples/project1_scar_cond1_H4K20me0_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,H4K20me0,project1_scar_cond1_strandedInput,1
project1_scar_cond1_strandedInput,/path/to/samples/project1_scar_cond1_strandedInput_R1.fastq.gz,/path/to/samples/project1_scar_cond1_strandedInput_R3.fastq.gz,/path/to/samples/project1_scar_cond1_strandedInput_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,,,
project1_scar_cond2_H3K9me3,/path/to/samples/project1_scar_cond2_H3K9me3_R1.fastq.gz,/path/to/samples/project1_scar_cond2_H3K9me3_R3.fastq.gz,/path/to/samples/project1_scar_cond2_H3K9me3_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,H3K9me3,project1_scar_cond2_strandedInput,1
project1_scar_cond2_H4K20me0,/path/to/samples/project1_scar_cond2_H4K20me0_R1.fastq.gz,/path/to/samples/project1_scar_cond2_H4K20me0_R3.fastq.gz,/path/to/samples/project1_scar_cond2_H4K20me0_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,H4K20me0,project1_scar_cond2_strandedInput,1
project1_scar_cond2_strandedInput,/path/to/samples/project1_scar_cond2_strandedInput_R1.fastq.gz,/path/to/samples/project1_scar_cond2_strandedInput_R3.fastq.gz,/path/to/samples/project1_scar_cond2_strandedInput_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,,,
project1_scar_cond3_H3K9me3,/path/to/samples/project1_scar_cond3_H3K9me3_R1.fastq.gz,/path/to/samples/project1_scar_cond3_H3K9me3_R3.fastq.gz,/path/to/samples/project1_scar_cond3_H3K9me3_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,H3K9me3,project1_scar_cond3_strandedInput,1
project1_scar_cond3_H4K20me0,/path/to/samples/project1_scar_cond3_H4K20me0_R1.fastq.gz,/path/to/samples/project1_scar_cond3_H4K20me0_R3.fastq.gz,/path/to/samples/project1_scar_cond3_H4K20me0_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,H4K20me0,project1_scar_cond3_strandedInput,1
project1_scar_cond3_strandedInput,/path/to/samples/project1_scar_cond3_strandedInput_R1.fastq.gz,/path/to/samples/project1_scar_cond3_strandedInput_R3.fastq.gz,/path/to/samples/project1_scar_cond3_strandedInput_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,,,
project1_scar_cond4_H3K9me3,/path/to/samples/project1_scar_cond4_H3K9me3_R1.fastq.gz,/path/to/samples/project1_scar_cond4_H3K9me3_R3.fastq.gz,/path/to/samples/project1_scar_cond4_H3K9me3_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,H3K9me3,project1_scar_cond4_strandedInput,1
project1_scar_cond4_strandedInput,/path/to/samples/project1_scar_cond4_strandedInput_R1.fastq.gz,/path/to/samples/project1_scar_cond4_strandedInput_R3.fastq.gz,/path/to/samples/project1_scar_cond4_strandedInput_R2.fastq.gz,/path/to/reference/OKseq_RFD_mESC_SRR7535256_R1.csorted.nodup.GRCm38_SE_smooth_results_w1000_s30_d30_z1.bed.gz,1,scarseq,reverse,,,
project1_chip_cond1_H3K9me3,/path/to/samples/project1_chip_cond1_H3K9me3_R1.fastq.gz,/path/to/samples/project1_chip_cond1_H3K9me3_R3.fastq.gz,/path/to/samples/project1_chip_cond1_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,project1_chip_cond1_Input,1
project1_chip_cond1_H3K9me3,/path/to/samples/project1_chip_cond1_H3K9me3_r2_R1.fastq.gz,/path/to/samples/project1_chip_cond1_H3K9me3_r2_R3.fastq.gz,/path/to/samples/project1_chip_cond1_H3K9me3_r2_R2.fastq.gz,,2,chipseq,,H3K9me3,project1_chip_cond1_Input,2
project1_chip_cond1_Input,/path/to/samples/project1_chip_cond1_Input_R1.fastq.gz,/path/to/samples/project1_chip_cond1_Input_R3.fastq.gz,/path/to/samples/project1_chip_cond1_Input_R2.fastq.gz,,1,chipseq,,,,
project1_chip_cond1_Input,/path/to/samples/project1_chip_cond1_Input_r2_R1.fastq.gz,/path/to/samples/project1_chip_cond1_Input_r2_R3.fastq.gz,/path/to/samples/project1_chip_cond1_Input_r2_R2.fastq.gz,,2,chipseq,,,,
project1_chip_cond2_H3K9me3,/path/to/samples/project1_chip_cond2_H3K9me3_R1.fastq.gz,/path/to/samples/project1_chip_cond2_H3K9me3_R3.fastq.gz,/path/to/samples/project1_chip_cond2_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,project1_chip_cond2_Input,1
project1_chip_cond2_H3K9me3,/path/to/samples/project1_chip_cond2_H3K9me3_r2_R1.fastq.gz,/path/to/samples/project1_chip_cond2_H3K9me3_r2_R3.fastq.gz,/path/to/samples/project1_chip_cond2_H3K9me3_r2_R2.fastq.gz,,2,chipseq,,H3K9me3,project1_chip_cond2_Input,2
project1_chip_cond2_Input,/path/to/samples/project1_chip_cond2_Input_R1.fastq.gz,/path/to/samples/project1_chip_cond2_Input_R3.fastq.gz,/path/to/samples/project1_chip_cond2_Input_R2.fastq.gz,,1,chipseq,,,,
project1_chip_cond2_Input,/path/to/samples/project1_chip_cond2_Input_r2_R1.fastq.gz,/path/to/samples/project1_chip_cond2_Input_r2_R3.fastq.gz,/path/to/samples/project1_chip_cond2_Input_r2_R2.fastq.gz,,2,chipseq,,,,
project1_chip_cond3_H3K9me3,/path/to/samples/project1_chip_cond3_H3K9me3_R1.fastq.gz,/path/to/samples/project1_chip_cond3_H3K9me3_R3.fastq.gz,/path/to/samples/project1_chip_cond3_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,project1_chip_cond3_Input,1
project1_chip_cond3_H3K9me3,/path/to/samples/project1_chip_cond3_H3K9me3_r2_R1.fastq.gz,/path/to/samples/project1_chip_cond3_H3K9me3_r2_R3.fastq.gz,/path/to/samples/project1_chip_cond3_H3K9me3_r2_R2.fastq.gz,,2,chipseq,,H3K9me3,project1_chip_cond3_Input,2
project1_chip_cond3_Input,/path/to/samples/project1_chip_cond3_Input_R1.fastq.gz,/path/to/samples/project1_chip_cond3_Input_R3.fastq.gz,/path/to/samples/project1_chip_cond3_Input_R2.fastq.gz,,1,chipseq,,,,
project1_chip_cond3_Input,/path/to/samples/project1_chip_cond3_Input_r2_R1.fastq.gz,/path/to/samples/project1_chip_cond3_Input_r2_R3.fastq.gz,/path/to/samples/project1_chip_cond3_Input_r2_R2.fastq.gz,,2,chipseq,,,,
project1_chip_cond4_H3K9me3,/path/to/samples/project1_chip_cond4_H3K9me3_R1.fastq.gz,/path/to/samples/project1_chip_cond4_H3K9me3_R3.fastq.gz,/path/to/samples/project1_chip_cond4_H3K9me3_R2.fastq.gz,,1,chipseq,,H3K9me3,project1_chip_cond4_Input,1
project1_chip_cond4_H3K9me3,/path/to/samples/project1_chip_cond4_H3K9me3_r2_R1.fastq.gz,/path/to/samples/project1_chip_cond4_H3K9me3_r2_R3.fastq.gz,/path/to/samples/project1_chip_cond4_H3K9me3_r2_R2.fastq.gz,,2,chipseq,,H3K9me3,project1_chip_cond4_Input,2
project1_chip_cond4_Input,/path/to/samples/project1_chip_cond4_Input_R1.fastq.gz,/path/to/samples/project1_chip_cond4_Input_R3.fastq.gz,/path/to/samples/project1_chip_cond4_Input_R2.fastq.gz,,1,chipseq,,,,
project1_chip_cond4_Input,/path/to/samples/project1_chip_cond4_Input_r2_R1.fastq.gz,/path/to/samples/project1_chip_cond4_Input_r2_R3.fastq.gz,/path/to/samples/project1_chip_cond4_Input_r2_R2.fastq.gz,,2,chipseq,,,,
```
</details>

### Reference genome files

#### Genome index, FASTA, GTF/GFF and gene BED files

The minimum reference genome requirements are a FASTA and a GTF file, all other files required to run the pipeline can be generated from these files. However, it is more storage and compute friendly if you are able to re-use reference genome files as efficiently as possible. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices (e.g. those unavailable on [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space. You can then either provide the appropriate reference genome files on the command-line via the appropriate parameters (e.g. `--bwa_index '/path/to/bwa/index/'`) or via a [custom config file](https://nf-co.re/usage/configuration#custom-configuration-files).

- If `--genome` is provided then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`.
- If `--gene_bed` is not provided then it will be generated from the GTF file.

#### Genome blacklist regions


#### OK-seq partitions

#### Initiation zones

## Credits

The [glseq](https://github.com/grothlab/glseq) pipeline was written by Samuel Ruiz-Pérez ([@samuelruizperez](https://github.com/samuelruizperez)) at the Groth Lab ([@grothlab](https://github.com/grothlab)).

Several scripts in this pipeline are based on [nf-core/chipseq](https://github.com/nf-core/chipseq) scripts, which were originally written by Chuan Wang ([@chuan-wang](https://github.com/chuan-wang)) and Phil Ewels ([@ewels](https://github.com/ewels)), re-implemented by Harshil Patel ([@drpatelh](https://github.com/drpatelh)), and converted to Nextflow DSL2 by Jose Espinosa-Carrasco ([@JoseEspinosa](https://github.com/JoseEspinosa)). For more information regarding nf-core/chipseq, see [https://github.com/nf-core/chipseq?tab=readme-ov-file#credits](https://github.com/nf-core/chipseq?tab=readme-ov-file#credits)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

<!-- TODO: -->
For further information or help, don't hesitate to get in touch through #######
## Citations

If you use [grothlab/glseq](https://github.com/grothlab/glseq) for your analysis, please cite it using the following DOI: #########

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

In addition, an extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
