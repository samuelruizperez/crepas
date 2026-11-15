# ChIP-seq

For this example, we will analyze the H3K27me3 and H3K9me3 ChIP-seq data of the same article as the [SCAR-seq tutorial](scarseq.md), in mouse embryonic stem cells grown in serum/LIF:

> Wenger, A., Biran, A., Alcaraz, N., et al. (2023). Symmetric inheritance of parental histones governs epigenome maintenance and embryonic stem cell identity. _Nature Genetics_, _55_, 1567–1578. [https://doi.org/10.1038/s41588-023-01476-x](https://doi.org/10.1038/s41588-023-01476-x)

These libraries have a _Drosophila melanogaster_ (dm6) spike-in, three biological replicates, and two technical replicates (sequencing runs) per biological replicate.

We start by creating a directory for this tutorial:

```bash
mkdir -p ${HOME}/projects/chipseq_tutorial
```

## Downloading the data

```bash
mkdir -p ${HOME}/projects/chipseq_tutorial/fetchngs/
```

### Input IDs file

The data is in GEO series [GSE154379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154379). Since each GEO sample has two sequencing runs, which we want to keep as technical replicates, we list the 72 SRA runs (SRR) instead of the GEO samples:

<details>
<summary>Show the 72 accessions</summary>

```csv
SRR23304630
SRR23304631
SRR23304632
SRR23304633
SRR23304634
SRR23304635
SRR23304636
SRR23304637
SRR23304638
SRR23304639
SRR23304640
SRR23304641
SRR23304642
SRR23304643
SRR23304644
SRR23304645
SRR23304646
SRR23304647
SRR23304648
SRR23304649
SRR23304650
SRR23304651
SRR23304652
SRR23304653
SRR23304654
SRR23304655
SRR23304656
SRR23304657
SRR23304658
SRR23304659
SRR23304660
SRR23304661
SRR23304662
SRR23304663
SRR23304664
SRR23304665
SRR23304666
SRR23304667
SRR23304668
SRR23304669
SRR23304670
SRR23304671
SRR23304672
SRR23304673
SRR23304674
SRR23304675
SRR23304676
SRR23304677
SRR23304678
SRR23304679
SRR23304680
SRR23304681
SRR23304682
SRR23304683
SRR23304684
SRR23304685
SRR23304686
SRR23304687
SRR23304688
SRR23304689
SRR23304690
SRR23304691
SRR23304692
SRR23304693
SRR23304694
SRR23304695
SRR23304696
SRR23304697
SRR23304698
SRR23304699
SRR23304700
SRR23304701
```

</details>

And we save it as `chipseq_tutorial.fetchngs.input.csv` in the `fetchngs` subdirectory.

### SBATCH script

We reuse the SBATCH script of the [OK-seq tutorial](okseq.md), replacing `okseq_tutorial` with `chipseq_tutorial` in the paths:

```bash
nextflow run nf-core/fetchngs \
    -r 1.12.0 \
    -profile singularity \
    --input ${HOME}/projects/chipseq_tutorial/fetchngs/chipseq_tutorial.fetchngs.input.csv \
    --outdir ${HOME}/projects/chipseq_tutorial/fetchngs/output/
```

and save it as `chipseq_tutorial.fetchngs.sbatch.sh`.

### Submitting the job

```bash
cd ${HOME}/projects/chipseq_tutorial/fetchngs/
sbatch chipseq_tutorial.fetchngs.sbatch.sh
```

## Running crepas

```bash
mkdir -p ${HOME}/projects/chipseq_tutorial/crepas/
```

### Input samplesheet CSV

Each row is one sequencing run. The two runs of the same library share `sample` and `biological_replicate`, and differ in `technical_replicate`; _crepas_ merges them after alignment. Each ChIP run points to the input run with the same biological and technical replicate (`input_control`, `input_control_biological_replicate` and `input_control_technical_replicate`). We save it as `chipseq_tutorial.crepas.samplesheet.csv`:

<details>
<summary>Show the samplesheet (72 rows)</summary>

```csv
sample,biological_replicate,technical_replicate,input_control,input_control_biological_replicate,input_control_technical_replicate,exp_type,antibody,fastq_1,fastq_2
ChIP_410_H3K27me3,1,1,Input_410_H3K27me3,1,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247442_SRR23304630_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247442_SRR23304630_2.fastq.gz
ChIP_410_H3K27me3,1,2,Input_410_H3K27me3,1,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247442_SRR23304631_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247442_SRR23304631_2.fastq.gz
ChIP_410_H3K27me3,2,1,Input_410_H3K27me3,2,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247445_SRR23304636_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247445_SRR23304636_2.fastq.gz
ChIP_410_H3K27me3,2,2,Input_410_H3K27me3,2,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247445_SRR23304637_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247445_SRR23304637_2.fastq.gz
ChIP_410_H3K27me3,3,1,Input_410_H3K27me3,3,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247448_SRR23304642_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247448_SRR23304642_2.fastq.gz
ChIP_410_H3K27me3,3,2,Input_410_H3K27me3,3,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247448_SRR23304643_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247448_SRR23304643_2.fastq.gz
Input_410_H3K27me3,1,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247451_SRR23304648_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247451_SRR23304648_2.fastq.gz
Input_410_H3K27me3,1,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247451_SRR23304649_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247451_SRR23304649_2.fastq.gz
Input_410_H3K27me3,2,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247454_SRR23304654_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247454_SRR23304654_2.fastq.gz
Input_410_H3K27me3,2,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247454_SRR23304655_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247454_SRR23304655_2.fastq.gz
Input_410_H3K27me3,3,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247457_SRR23304660_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247457_SRR23304660_2.fastq.gz
Input_410_H3K27me3,3,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247457_SRR23304661_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247457_SRR23304661_2.fastq.gz
ChIP_442_H3K27me3,1,1,Input_442_H3K27me3,1,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247443_SRR23304632_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247443_SRR23304632_2.fastq.gz
ChIP_442_H3K27me3,1,2,Input_442_H3K27me3,1,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247443_SRR23304633_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247443_SRR23304633_2.fastq.gz
ChIP_442_H3K27me3,2,1,Input_442_H3K27me3,2,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247446_SRR23304638_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247446_SRR23304638_2.fastq.gz
ChIP_442_H3K27me3,2,2,Input_442_H3K27me3,2,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247446_SRR23304639_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247446_SRR23304639_2.fastq.gz
ChIP_442_H3K27me3,3,1,Input_442_H3K27me3,3,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247449_SRR23304644_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247449_SRR23304644_2.fastq.gz
ChIP_442_H3K27me3,3,2,Input_442_H3K27me3,3,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247449_SRR23304645_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247449_SRR23304645_2.fastq.gz
Input_442_H3K27me3,1,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247452_SRR23304650_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247452_SRR23304650_2.fastq.gz
Input_442_H3K27me3,1,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247452_SRR23304651_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247452_SRR23304651_2.fastq.gz
Input_442_H3K27me3,2,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247455_SRR23304656_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247455_SRR23304656_2.fastq.gz
Input_442_H3K27me3,2,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247455_SRR23304657_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247455_SRR23304657_2.fastq.gz
Input_442_H3K27me3,3,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247458_SRR23304662_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247458_SRR23304662_2.fastq.gz
Input_442_H3K27me3,3,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247458_SRR23304663_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247458_SRR23304663_2.fastq.gz
ChIP_588_H3K27me3,1,1,Input_588_H3K27me3,1,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247444_SRR23304634_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247444_SRR23304634_2.fastq.gz
ChIP_588_H3K27me3,1,2,Input_588_H3K27me3,1,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247444_SRR23304635_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247444_SRR23304635_2.fastq.gz
ChIP_588_H3K27me3,2,1,Input_588_H3K27me3,2,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247447_SRR23304640_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247447_SRR23304640_2.fastq.gz
ChIP_588_H3K27me3,2,2,Input_588_H3K27me3,2,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247447_SRR23304641_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247447_SRR23304641_2.fastq.gz
ChIP_588_H3K27me3,3,1,Input_588_H3K27me3,3,1,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247450_SRR23304646_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247450_SRR23304646_2.fastq.gz
ChIP_588_H3K27me3,3,2,Input_588_H3K27me3,3,2,ChIP-seq,H3K27me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247450_SRR23304647_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247450_SRR23304647_2.fastq.gz
Input_588_H3K27me3,1,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247453_SRR23304652_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247453_SRR23304652_2.fastq.gz
Input_588_H3K27me3,1,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247453_SRR23304653_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247453_SRR23304653_2.fastq.gz
Input_588_H3K27me3,2,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247456_SRR23304658_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247456_SRR23304658_2.fastq.gz
Input_588_H3K27me3,2,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247456_SRR23304659_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247456_SRR23304659_2.fastq.gz
Input_588_H3K27me3,3,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247459_SRR23304664_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247459_SRR23304664_2.fastq.gz
Input_588_H3K27me3,3,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247459_SRR23304665_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247459_SRR23304665_2.fastq.gz
ChIP_410_H3K9me3,1,1,Input_410_H3K9me3,1,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247460_SRR23304666_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247460_SRR23304666_2.fastq.gz
ChIP_410_H3K9me3,1,2,Input_410_H3K9me3,1,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247460_SRR23304667_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247460_SRR23304667_2.fastq.gz
ChIP_410_H3K9me3,2,1,Input_410_H3K9me3,2,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247463_SRR23304672_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247463_SRR23304672_2.fastq.gz
ChIP_410_H3K9me3,2,2,Input_410_H3K9me3,2,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247463_SRR23304673_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247463_SRR23304673_2.fastq.gz
ChIP_410_H3K9me3,3,1,Input_410_H3K9me3,3,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247466_SRR23304678_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247466_SRR23304678_2.fastq.gz
ChIP_410_H3K9me3,3,2,Input_410_H3K9me3,3,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247466_SRR23304679_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247466_SRR23304679_2.fastq.gz
Input_410_H3K9me3,1,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247469_SRR23304684_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247469_SRR23304684_2.fastq.gz
Input_410_H3K9me3,1,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247469_SRR23304685_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247469_SRR23304685_2.fastq.gz
Input_410_H3K9me3,2,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247472_SRR23304690_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247472_SRR23304690_2.fastq.gz
Input_410_H3K9me3,2,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247472_SRR23304691_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247472_SRR23304691_2.fastq.gz
Input_410_H3K9me3,3,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247475_SRR23304696_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247475_SRR23304696_2.fastq.gz
Input_410_H3K9me3,3,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247475_SRR23304697_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247475_SRR23304697_2.fastq.gz
ChIP_442_H3K9me3,1,1,Input_442_H3K9me3,1,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247461_SRR23304668_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247461_SRR23304668_2.fastq.gz
ChIP_442_H3K9me3,1,2,Input_442_H3K9me3,1,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247461_SRR23304669_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247461_SRR23304669_2.fastq.gz
ChIP_442_H3K9me3,2,1,Input_442_H3K9me3,2,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247464_SRR23304674_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247464_SRR23304674_2.fastq.gz
ChIP_442_H3K9me3,2,2,Input_442_H3K9me3,2,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247464_SRR23304675_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247464_SRR23304675_2.fastq.gz
ChIP_442_H3K9me3,3,1,Input_442_H3K9me3,3,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247467_SRR23304680_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247467_SRR23304680_2.fastq.gz
ChIP_442_H3K9me3,3,2,Input_442_H3K9me3,3,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247467_SRR23304681_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247467_SRR23304681_2.fastq.gz
Input_442_H3K9me3,1,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247470_SRR23304686_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247470_SRR23304686_2.fastq.gz
Input_442_H3K9me3,1,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247470_SRR23304687_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247470_SRR23304687_2.fastq.gz
Input_442_H3K9me3,2,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247473_SRR23304692_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247473_SRR23304692_2.fastq.gz
Input_442_H3K9me3,2,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247473_SRR23304693_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247473_SRR23304693_2.fastq.gz
Input_442_H3K9me3,3,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247476_SRR23304698_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247476_SRR23304698_2.fastq.gz
Input_442_H3K9me3,3,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247476_SRR23304699_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247476_SRR23304699_2.fastq.gz
ChIP_588_H3K9me3,1,1,Input_588_H3K9me3,1,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247462_SRR23304670_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247462_SRR23304670_2.fastq.gz
ChIP_588_H3K9me3,1,2,Input_588_H3K9me3,1,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247462_SRR23304671_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247462_SRR23304671_2.fastq.gz
ChIP_588_H3K9me3,2,1,Input_588_H3K9me3,2,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247465_SRR23304676_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247465_SRR23304676_2.fastq.gz
ChIP_588_H3K9me3,2,2,Input_588_H3K9me3,2,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247465_SRR23304677_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247465_SRR23304677_2.fastq.gz
ChIP_588_H3K9me3,3,1,Input_588_H3K9me3,3,1,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247468_SRR23304682_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247468_SRR23304682_2.fastq.gz
ChIP_588_H3K9me3,3,2,Input_588_H3K9me3,3,2,ChIP-seq,H3K9me3,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247468_SRR23304683_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247468_SRR23304683_2.fastq.gz
Input_588_H3K9me3,1,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247471_SRR23304688_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247471_SRR23304688_2.fastq.gz
Input_588_H3K9me3,1,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247471_SRR23304689_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247471_SRR23304689_2.fastq.gz
Input_588_H3K9me3,2,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247474_SRR23304694_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247474_SRR23304694_2.fastq.gz
Input_588_H3K9me3,2,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247474_SRR23304695_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247474_SRR23304695_2.fastq.gz
Input_588_H3K9me3,3,1,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247477_SRR23304700_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247477_SRR23304700_2.fastq.gz
Input_588_H3K9me3,3,2,,,,ChIP-seq,,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247477_SRR23304701_1.fastq.gz,<PATH>/projects/chipseq_tutorial/fetchngs/output/fastq/SRX19247477_SRR23304701_2.fastq.gz
```

</details>

> [!NOTE]
> The sample names follow the GEO sample names: `410`, `442` and `588` are the three cell lines of the article, and the replicate numbers are those on GEO.

### Parameter YAML file

For spike-in experiments, _crepas_ aligns the reads to a hybrid genome with both the endogenous (mouse) and the spike-in (fly) chromosomes, and `spikein_genome` gives the suffix of the spike-in chromosome names (e.g. `chr2L_dm6`). We provide this hybrid genome with `hybrid_fasta`, instead of `fasta`. We also enable TE counting, as in the [SCAR-seq tutorial](scarseq.md):

```yaml
input: <PATH>/projects/chipseq_tutorial/crepas/chipseq_tutorial.crepas.samplesheet.csv
outdir: <PATH>/projects/chipseq_tutorial/crepas/output
multiqc_title: ChIP-seq_H3K27me3_H3K9me3
genome: GRCm39
aligner: chromap
hybrid_fasta: <GRCm39_DM6_HYBRID_FASTA>
spikein_genome: dm6
gtf: <GRCm39_GTF>
blacklist: <GRCm39_BLACKLIST_BED>
te_gtf: <GRCm39_TE_GTF>
seq_platform: ILLUMINA
read_length: 100
fragment_size: 1000
coverage_extend_reads: 1000
coverage_bin_size: 1000
map_n_multimappers: 100
save_align_intermeds: true
genes_plotprofile_mode: reference_point
skip_plot_fingerprint: true
skip_spp_qc: true
skip_flT3: false
skip_te_counting: false
skip_telocal: false
```

Save this file as `chipseq_tutorial.crepas.params.yaml` in the `crepas/` subdirectory.

### SBATCH script

```bash
cd ${HOME}/projects/chipseq_tutorial/crepas/

nextflow run grothlab/crepas \
    -r 1.2.0 \
    -profile singularity \
    -params-file ${HOME}/projects/chipseq_tutorial/crepas/chipseq_tutorial.crepas.params.yaml
```

Save it as `chipseq_tutorial.crepas.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/chipseq_tutorial/crepas/
sbatch chipseq_tutorial.crepas.sbatch.sh
```
