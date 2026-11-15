# ChOR-seq

For this example, we will analyze the H2BK120ub ChOR-seq time course from the following article:

> Flury, V., Reverón-Gómez, N., Alcaraz, N., Stewart-Morgan, K. R., Wenger, A., Klose, R. J., & Groth, A. (2023). Recycling of modified H2A-H2B provides short-term memory of chromatin states. _Cell_, _186_(5), 1050–1065.e19. [https://doi.org/10.1016/j.cell.2023.01.007](https://doi.org/10.1016/j.cell.2023.01.007)

The ChOR-seq libraries were made at several times after the nascent DNA labelling, each with its own clicked input, and also have a _Drosophila melanogaster_ (dm6) spike-in.

We start by creating a directory for this tutorial:

```bash
mkdir -p ${HOME}/projects/chorseq_tutorial
```

## Downloading the data

```bash
mkdir -p ${HOME}/projects/chorseq_tutorial/fetchngs/
```

### Input IDs file

The data is in GEO series [GSE204981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204981). We list the GEO samples of H2BK120ub and their clicked inputs:

<details>
<summary>Show the 26 accessions</summary>

```csv
GSM6203514
GSM6203515
GSM6203516
GSM6203517
GSM6203518
GSM6203519
GSM6203520
GSM6203521
GSM6203522
GSM6203523
GSM6203524
GSM6203525
GSM6203526
GSM6203527
GSM6203534
GSM6203535
GSM6203541
GSM6203544
GSM6203546
GSM6203547
GSM6203548
GSM6203549
GSM6203550
GSM6203551
GSM6203538
GSM6203539
```

</details>

And we save it as `chorseq_tutorial.fetchngs.input.csv` in the `fetchngs` subdirectory.

### SBATCH script

```bash
nextflow run nf-core/fetchngs \
    -r 1.12.0 \
    -profile singularity \
    --input ${HOME}/projects/chorseq_tutorial/fetchngs/chorseq_tutorial.fetchngs.input.csv \
    --outdir ${HOME}/projects/chorseq_tutorial/fetchngs/output/
```

Save it as `chorseq_tutorial.fetchngs.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/chorseq_tutorial/fetchngs/
sbatch chorseq_tutorial.fetchngs.sbatch.sh
```

## Running crepas

```bash
mkdir -p ${HOME}/projects/chorseq_tutorial/crepas/
```

### Input samplesheet CSV

`T0` to `T480` in the sample names are the minutes after the labelling. Not every time point has a clicked input for both replicates, so some ChOR-seq samples use the input of replicate 1 for both of their replicates. We save it as `chorseq_tutorial.crepas.samplesheet.csv`:

<details>
<summary>Show the samplesheet (26 rows)</summary>

```csv
sample,biological_replicate,input_control,input_control_biological_replicate,exp_type,antibody,fastq_1
ChOR_T0_H2BK120ub,1,ChOR_T0_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725677_SRR19676829.fastq.gz
ChOR_T0_H2BK120ub,2,ChOR_T0_H2BK120ub_clickedInput,2,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725678_SRR19676828.fastq.gz
ChOR_T15_H2BK120ub,1,ChOR_T15_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725679_SRR19676827.fastq.gz
ChOR_T15_H2BK120ub,2,ChOR_T15_H2BK120ub_clickedInput,2,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725680_SRR19676826.fastq.gz
ChOR_T30_H2BK120ub,1,ChOR_T30_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725681_SRR19676825.fastq.gz
ChOR_T30_H2BK120ub,2,ChOR_T30_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725682_SRR19676824.fastq.gz
ChOR_T60_H2BK120ub,1,ChOR_T60_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725683_SRR19676823.fastq.gz
ChOR_T60_H2BK120ub,2,ChOR_T60_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725684_SRR19676822.fastq.gz
ChOR_T120_H2BK120ub,1,ChOR_T120_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725685_SRR19676821.fastq.gz
ChOR_T120_H2BK120ub,2,ChOR_T120_H2BK120ub_clickedInput,2,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725686_SRR19676820.fastq.gz
ChOR_T180_H2BK120ub,1,ChOR_T180_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725687_SRR19676819.fastq.gz
ChOR_T180_H2BK120ub,2,ChOR_T180_H2BK120ub_clickedInput,2,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725688_SRR19676818.fastq.gz
ChOR_T480_H2BK120ub,1,ChOR_T480_H2BK120ub_clickedInput,1,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725689_SRR19676817.fastq.gz
ChOR_T480_H2BK120ub,2,ChOR_T480_H2BK120ub_clickedInput,2,ChOR-seq,H2BK120ub,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725690_SRR19676816.fastq.gz
ChOR_T0_H2BK120ub_clickedInput,1,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725697_SRR19676809.fastq.gz
ChOR_T0_H2BK120ub_clickedInput,2,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725698_SRR19676808.fastq.gz
ChOR_T15_H2BK120ub_clickedInput,1,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725701_SRR19676805.fastq.gz
ChOR_T15_H2BK120ub_clickedInput,2,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725702_SRR19676804.fastq.gz
ChOR_T30_H2BK120ub_clickedInput,1,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725704_SRR19676802.fastq.gz
ChOR_T60_H2BK120ub_clickedInput,1,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725723_SRR19676783.fastq.gz
ChOR_T120_H2BK120ub_clickedInput,1,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725725_SRR19676781.fastq.gz
ChOR_T120_H2BK120ub_clickedInput,2,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725726_SRR19676780.fastq.gz
ChOR_T180_H2BK120ub_clickedInput,1,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725727_SRR19676779.fastq.gz
ChOR_T180_H2BK120ub_clickedInput,2,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725728_SRR19676778.fastq.gz
ChOR_T480_H2BK120ub_clickedInput,1,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725729_SRR19676777.fastq.gz
ChOR_T480_H2BK120ub_clickedInput,2,,,ChOR-seq,,<PATH>/projects/chorseq_tutorial/fetchngs/output/fastq/SRX15725730_SRR19676776.fastq.gz
```

</details>

### Parameter YAML file

As in the [ChIP-seq tutorial](chipseq.md), we give the hybrid mouse and fly genome with `hybrid_fasta` and enable TE counting:

```yaml
input: <PATH>/projects/chorseq_tutorial/crepas/chorseq_tutorial.crepas.samplesheet.csv
outdir: <PATH>/projects/chorseq_tutorial/crepas/output
multiqc_title: ChOR-seq_H2BK120ub
genome: GRCm39
aligner: chromap
hybrid_fasta: <GRCm39_DM6_HYBRID_FASTA>
spikein_genome: dm6
gtf: <GRCm39_GTF>
blacklist: <GRCm39_BLACKLIST_BED>
te_gtf: <GRCm39_TE_GTF>
seq_platform: ILLUMINA
read_length: 75
fragment_size: 250
coverage_extend_reads: 250
coverage_bin_size: 50
map_n_multimappers: 100
save_align_intermeds: true
skip_plot_fingerprint: true
skip_spp_qc: true
skip_flT3: true
skip_te_counting: false
skip_telocal: false
```

Save this file as `chorseq_tutorial.crepas.params.yaml` in the `crepas/` subdirectory.

### SBATCH script

```bash
cd ${HOME}/projects/chorseq_tutorial/crepas/

nextflow run grothlab/crepas \
    -r 1.2.0 \
    -profile singularity \
    -params-file ${HOME}/projects/chorseq_tutorial/crepas/chorseq_tutorial.crepas.params.yaml
```

Save it as `chorseq_tutorial.crepas.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/chorseq_tutorial/crepas/
sbatch chorseq_tutorial.crepas.sbatch.sh
```
