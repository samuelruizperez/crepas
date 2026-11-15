# Introduction

These tutorials show how to analyze public datasets of different experiment types with _crepas_, from downloading the raw data to the pipeline results:

1. [OK-seq](okseq.md): replication fork directionality and initiation zones in mouse embryonic stem cells.
2. [SCAR-seq](scarseq.md): partitioning of H3K27me3 and H3K27ac between sister chromatids, using the results of the OK-seq tutorial.
3. [ChIP-seq](chipseq.md): H3K27me3 and H3K9me3 with a spike-in and technical replicates.
4. [ChOR-seq](chorseq.md): a time course of H2BK120ub on nascent chromatin, with a spike-in.
5. [eSPAN](espan.md): H3K4me3 partitioning in fission yeast mutants.
6. [Early/Late Repli-seq](el_repliseq.md): early/mid/late replication timing of embryonic stem cells and 2-cell-like cells.
7. [High-resolution Repli-seq](hr_repliseq.md): 16-fraction replication timing of human embryonic stem cells.

## Setting the environment

Before analyzing any kind of data, we need to install or load Nextflow and a container runtime, as explained in the [Quick start section of the README](../../../README.md#quick-start).

> [!NOTE]
> In this case, we will use [_SBATCH (SLURM)_](https://slurm.schedmd.com/documentation.html) to schedule jobs,
> [_modules_](https://modules.readthedocs.io/en/latest/) to load Nextflow, and [_singularity_](https://sylabs.io/guides/3.0/user-guide/)
> to run the containerized software in the pipeline. However, these choices
> may vary depending on your specific computing environment, so you should check your local documentation.

The parameter files of the tutorials refer to reference files with placeholders such as `<GRCm39_FASTA>` or `<GRCm39_BLACKLIST_BED>`. Replace them with the **full paths** to your own copies of these files (see [Reference genome files](../../usage.md#reference-genome-files)). _crepas_ builds the aligner index from the FASTA file in the first run and, with `--save_reference` (enabled by default), saves it in the output directory, so that later runs can reuse it with `--chromap_index`.
