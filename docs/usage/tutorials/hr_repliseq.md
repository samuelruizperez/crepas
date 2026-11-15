# High-resolution Repli-seq

For this example, we will analyze the high-resolution Repli-seq data of H1 human embryonic stem cells from the following article:

> Zhao, P. A., Sasaki, T., & Gilbert, D. M. (2020). High-resolution Repli-Seq defines the temporal choreography of initiation, elongation and termination of replication in mammalian cells. _Genome Biology_, _21_, 76. [https://doi.org/10.1186/s13059-020-01983-8](https://doi.org/10.1186/s13059-020-01983-8)

In high-resolution Repli-seq, S phase is sorted into 16 fractions (`S1` to `S16`), and a `G1` library is used as a control for mappability. In the same run, we will also analyze an early/late (E/L) Repli-seq dataset of H1 cells from the [4D Nucleome data portal](https://data.4dnucleome.org/), to compare the two designs.

We start by creating a directory for this tutorial:

```bash
mkdir -p ${HOME}/projects/hr_repliseq_tutorial
```

## Downloading the data

```bash
mkdir -p ${HOME}/projects/hr_repliseq_tutorial/fetchngs/
```

### Input IDs file

The high-resolution Repli-seq data of H1 cells are 17 runs of SRA project [PRJNA566530](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA566530):

<details>
<summary>Show the 17 accessions</summary>

```csv
SRR10150004
SRR10150005
SRR10150006
SRR10150007
SRR10150008
SRR10150009
SRR10150010
SRR10150011
SRR10150012
SRR10150013
SRR10150014
SRR10150015
SRR10150016
SRR10150017
SRR10150018
SRR10150019
SRR10150020
```

</details>

And we save it as `hr_repliseq_tutorial.fetchngs.input.csv` in the `fetchngs` subdirectory.

### SBATCH script

```bash
nextflow run nf-core/fetchngs \
    -r 1.12.0 \
    -profile singularity \
    --input ${HOME}/projects/hr_repliseq_tutorial/fetchngs/hr_repliseq_tutorial.fetchngs.input.csv \
    --outdir ${HOME}/projects/hr_repliseq_tutorial/fetchngs/output/
```

Save it as `hr_repliseq_tutorial.fetchngs.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/hr_repliseq_tutorial/fetchngs/
sbatch hr_repliseq_tutorial.fetchngs.sbatch.sh
```

### Downloading the E/L Repli-seq data

The E/L Repli-seq files are not on SRA, so we download them directly from the 4D Nucleome data portal:

```bash
mkdir -p ${HOME}/projects/hr_repliseq_tutorial/4dn/
cd ${HOME}/projects/hr_repliseq_tutorial/4dn/
wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/077e1e46-792c-43e7-8f5c-36316a129ae1/4DNFIVJP3AYA.fastq.gz
wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/4685aedb-29f1-4d58-b4a1-7ac2ee5ecc9c/4DNFIJLCP2D2.fastq.gz
wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/8f556d77-5992-4b31-ba97-985127b3e12d/4DNFIZUZHKUP.fastq.gz
wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/041f92dc-95fd-433c-94bf-ed938cdd0334/4DNFIGIS77SL.fastq.gz
wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/bc9de820-f133-42af-bbee-0d335e846e41/4DNFIPBQPYG5.fastq.gz
wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/d1d3531f-e06a-4b0a-aae1-d0130b965730/4DNFIPL7BXYA.fastq.gz
```

## Running crepas

```bash
mkdir -p ${HOME}/projects/hr_repliseq_tutorial/crepas/
```

### Input samplesheet CSV

The 16 S-phase fractions and the `G1` control of the high-resolution design, and the `early` and `late` fractions of the E/L design, all go in the same samplesheet with the same `sample` name; _crepas_ analyzes each design separately. In the E/L data, replicate 2 was sequenced twice, so it has two rows per fraction. We save it as `hr_repliseq_tutorial.crepas.samplesheet.csv`:

<details>
<summary>Show the samplesheet (23 rows)</summary>

```csv
sample,biological_replicate,exp_type,rt_fraction,fastq_1
H1,1,Repli-seq,G1,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875434_SRR10150004.fastq.gz
H1,1,Repli-seq,S1,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875435_SRR10150005.fastq.gz
H1,1,Repli-seq,S2,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875436_SRR10150006.fastq.gz
H1,1,Repli-seq,S3,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875437_SRR10150007.fastq.gz
H1,1,Repli-seq,S4,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875438_SRR10150008.fastq.gz
H1,1,Repli-seq,S5,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875439_SRR10150009.fastq.gz
H1,1,Repli-seq,S6,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875440_SRR10150010.fastq.gz
H1,1,Repli-seq,S7,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875441_SRR10150011.fastq.gz
H1,1,Repli-seq,S8,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875442_SRR10150012.fastq.gz
H1,1,Repli-seq,S9,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875443_SRR10150013.fastq.gz
H1,1,Repli-seq,S10,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875444_SRR10150014.fastq.gz
H1,1,Repli-seq,S11,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875445_SRR10150015.fastq.gz
H1,1,Repli-seq,S12,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875446_SRR10150016.fastq.gz
H1,1,Repli-seq,S13,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875447_SRR10150017.fastq.gz
H1,1,Repli-seq,S14,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875448_SRR10150018.fastq.gz
H1,1,Repli-seq,S15,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875449_SRR10150019.fastq.gz
H1,1,Repli-seq,S16,<PATH>/projects/hr_repliseq_tutorial/fetchngs/output/fastq/SRX6875450_SRR10150020.fastq.gz
H1,1,Repli-seq,early,<PATH>/projects/hr_repliseq_tutorial/4dn/4DNFIVJP3AYA.fastq.gz
H1,2,Repli-seq,early,<PATH>/projects/hr_repliseq_tutorial/4dn/4DNFIJLCP2D2.fastq.gz
H1,2,Repli-seq,early,<PATH>/projects/hr_repliseq_tutorial/4dn/4DNFIZUZHKUP.fastq.gz
H1,1,Repli-seq,late,<PATH>/projects/hr_repliseq_tutorial/4dn/4DNFIGIS77SL.fastq.gz
H1,2,Repli-seq,late,<PATH>/projects/hr_repliseq_tutorial/4dn/4DNFIPBQPYG5.fastq.gz
H1,2,Repli-seq,late,<PATH>/projects/hr_repliseq_tutorial/4dn/4DNFIPL7BXYA.fastq.gz
```

</details>

### Parameter YAML file

```yaml
input: <PATH>/projects/hr_repliseq_tutorial/crepas/hr_repliseq_tutorial.crepas.samplesheet.csv
outdir: <PATH>/projects/hr_repliseq_tutorial/crepas/output
multiqc_title: Repli-seq_H1
genome: GRCh38
aligner: chromap
fasta: <GRCh38_FASTA>
gtf: <GRCh38_GTF>
blacklist: <GRCh38_BLACKLIST_BED>
read_length: 50
skip_peak_calling: true
skip_flTbl: true
fingerprint_bins: 200
```

Save this file as `hr_repliseq_tutorial.crepas.params.yaml` in the `crepas/` subdirectory.

### SBATCH script

```bash
cd ${HOME}/projects/hr_repliseq_tutorial/crepas/

nextflow run grothlab/crepas \
    -r 1.2.0 \
    -profile singularity \
    -params-file ${HOME}/projects/hr_repliseq_tutorial/crepas/hr_repliseq_tutorial.crepas.params.yaml
```

Save it as `hr_repliseq_tutorial.crepas.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/hr_repliseq_tutorial/crepas/
sbatch hr_repliseq_tutorial.crepas.sbatch.sh
```

The high-resolution Repli-seq results are in `crepas/output/chromap/mergedLibrary/Repli-seq/hr_repliseq/`, and the E/L results in `crepas/output/chromap/mergedLibrary/Repli-seq/EL_repliseq/`.
