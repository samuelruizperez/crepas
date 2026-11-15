# Early/Late Repli-seq

For this example, we will compare the replication timing of mouse embryonic stem cells (ESCs) and 2-cell-like cells (2CLCs) from the following article:

> Nakatani, T., et al. (2022). DNA replication fork speed underlies cell fate changes and promotes reprogramming. _Nature Genetics_, _54_(3), 318–327. [https://doi.org/10.1038/s41588-022-01023-0](https://doi.org/10.1038/s41588-022-01023-0)

The cells were sorted into early, mid and late S-phase fractions, with two biological replicates each.

We start by creating a directory for this tutorial:

```bash
mkdir -p ${HOME}/projects/el_repliseq_tutorial
```

## Downloading the data

```bash
mkdir -p ${HOME}/projects/el_repliseq_tutorial/fetchngs/
```

### Input IDs file

The data is in SRA project [PRJNA561657](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA561657). We list its 12 SRA runs:

<details>
<summary>Show the 12 accessions</summary>

```csv
SRR10022299
SRR10022300
SRR10022301
SRR10022302
SRR10022303
SRR10022304
SRR10022305
SRR10022306
SRR10022307
SRR10022308
SRR10022309
SRR10022310
```

</details>

And we save it as `el_repliseq_tutorial.fetchngs.input.csv` in the `fetchngs` subdirectory.

### SBATCH script

```bash
nextflow run nf-core/fetchngs \
    -r 1.12.0 \
    -profile singularity \
    --input ${HOME}/projects/el_repliseq_tutorial/fetchngs/el_repliseq_tutorial.fetchngs.input.csv \
    --outdir ${HOME}/projects/el_repliseq_tutorial/fetchngs/output/
```

Save it as `el_repliseq_tutorial.fetchngs.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/el_repliseq_tutorial/fetchngs/
sbatch el_repliseq_tutorial.fetchngs.sbatch.sh
```

## Running crepas

```bash
mkdir -p ${HOME}/projects/el_repliseq_tutorial/crepas/
```

### Input samplesheet CSV

For Repli-seq, the `rt_fraction` column gives the S-phase fraction of each library (`early`, `mid` or `late`), and the libraries of the same cells share the `sample` name. On SRA, the ESCs are the `GFP-` and the 2CLCs the `GFP+` samples. We save it as `el_repliseq_tutorial.crepas.samplesheet.csv`:

<details>
<summary>Show the samplesheet (12 rows)</summary>

```csv
sample,biological_replicate,exp_type,rt_fraction,fastq_1
mESC,1,Repli-seq,early,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759360_SRR10022299.fastq.gz
mESC,2,Repli-seq,early,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759361_SRR10022300.fastq.gz
mESC,1,Repli-seq,mid,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759362_SRR10022301.fastq.gz
mESC,2,Repli-seq,mid,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759363_SRR10022302.fastq.gz
mESC,1,Repli-seq,late,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759364_SRR10022303.fastq.gz
mESC,2,Repli-seq,late,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759365_SRR10022304.fastq.gz
2CLCs,1,Repli-seq,early,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759366_SRR10022305.fastq.gz
2CLCs,2,Repli-seq,early,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759367_SRR10022306.fastq.gz
2CLCs,1,Repli-seq,mid,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759368_SRR10022307.fastq.gz
2CLCs,2,Repli-seq,mid,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759369_SRR10022308.fastq.gz
2CLCs,1,Repli-seq,late,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759370_SRR10022309.fastq.gz
2CLCs,2,Repli-seq,late,<PATH>/projects/el_repliseq_tutorial/fetchngs/output/fastq/SRX6759371_SRR10022310.fastq.gz
```

</details>

### Parameter YAML file

Repli-seq does not need peak calling, so we skip it:

```yaml
input: <PATH>/projects/el_repliseq_tutorial/crepas/el_repliseq_tutorial.crepas.samplesheet.csv
outdir: <PATH>/projects/el_repliseq_tutorial/crepas/output
multiqc_title: Repli-seq_ESC_2CLC
genome: GRCm39
aligner: chromap
fasta: <GRCm39_FASTA>
gtf: <GRCm39_GTF>
blacklist: <GRCm39_BLACKLIST_BED>
read_length: 75
skip_peak_calling: true
skip_flTbl: true
fingerprint_bins: 200
```

Save this file as `el_repliseq_tutorial.crepas.params.yaml` in the `crepas/` subdirectory.

### SBATCH script

```bash
cd ${HOME}/projects/el_repliseq_tutorial/crepas/

nextflow run grothlab/crepas \
    -r 1.2.0 \
    -profile singularity \
    -params-file ${HOME}/projects/el_repliseq_tutorial/crepas/el_repliseq_tutorial.crepas.params.yaml
```

Save it as `el_repliseq_tutorial.crepas.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/el_repliseq_tutorial/crepas/
sbatch el_repliseq_tutorial.crepas.sbatch.sh
```

The replication timing tracks are in `crepas/output/chromap/mergedLibrary/Repli-seq/EL_repliseq/`. See the [output documentation](../../output.md) for the description of each of them.
