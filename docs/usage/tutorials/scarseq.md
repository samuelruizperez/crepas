# SCAR-seq

For this example, we will replicate the SCAR-seq analyses of H3K27me3 and H3K27ac in wild-type and MCM2-2A mouse embryonic stem cells from the following article:

> Wenger, A., Biran, A., Alcaraz, N., et al. (2023). Symmetric inheritance of parental histones governs epigenome maintenance and embryonic stem cell identity. _Nature Genetics_, _55_, 1567–1578. [https://doi.org/10.1038/s41588-023-01476-x](https://doi.org/10.1038/s41588-023-01476-x)

SCAR-seq partitions are plotted around replication initiation zones, and _crepas_ takes these zones and the replication fork directionality from an OK-seq run. Therefore, we will use the outputs of the [OK-seq tutorial](okseq.md), so run that one first.

We start by creating a directory for this tutorial:

```bash
mkdir -p ${HOME}/projects/scarseq_tutorial
```

## Downloading the data

As in the [OK-seq tutorial](okseq.md), we download the data with nf-core/fetchngs:

```bash
mkdir -p ${HOME}/projects/scarseq_tutorial/fetchngs/
```

### Input IDs file

The SCAR-seq libraries and their stranded inputs are in GEO series [GSE154380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154380). nf-core/fetchngs also accepts GEO sample IDs (GSM), so we list the 56 samples we need:

<details>
<summary>Show the 56 accessions</summary>

```csv
GSM4670097
GSM4670098
GSM4670099
GSM4670103
GSM4670104
GSM4670105
GSM4670106
GSM4670107
GSM4670108
GSM4670112
GSM4670113
GSM4670114
GSM4670133
GSM4670134
GSM4670135
GSM4670139
GSM4670140
GSM4670141
GSM5858868
GSM5858869
GSM5858870
GSM5858871
GSM5858872
GSM5858873
GSM5858874
GSM5858875
GSM5858880
GSM5858881
GSM5858882
GSM5858883
GSM5858884
GSM5858885
GSM5858886
GSM5858887
GSM5858892
GSM5858893
GSM5858894
GSM5858895
GSM5858896
GSM5858897
GSM5858898
GSM5858899
GSM5858904
GSM5858905
GSM5858906
GSM5858907
GSM5858908
GSM5858909
GSM5858910
GSM5858911
GSM5858925
GSM5858926
GSM5858927
GSM5858931
GSM5858932
GSM5858933
```

</details>

And we save it as `scarseq_tutorial.fetchngs.input.csv` in the `fetchngs` subdirectory.

### SBATCH script

We can reuse the SBATCH script of the [OK-seq tutorial](okseq.md), replacing `okseq_tutorial` with `scarseq_tutorial` in the paths:

```bash
nextflow run nf-core/fetchngs \
    -r 1.12.0 \
    -profile singularity \
    --input ${HOME}/projects/scarseq_tutorial/fetchngs/scarseq_tutorial.fetchngs.input.csv \
    --outdir ${HOME}/projects/scarseq_tutorial/fetchngs/output/
```

and save it as `scarseq_tutorial.fetchngs.sbatch.sh`.

### Submitting the job

```bash
cd ${HOME}/projects/scarseq_tutorial/fetchngs/
sbatch scarseq_tutorial.fetchngs.sbatch.sh
```

When this job is finished, the FastQ files are in the `fetchngs/output/fastq/` subdirectory, named `<SRX>_<SRR>.fastq.gz` for single-end and `<SRX>_<SRR>_1.fastq.gz`/`<SRX>_<SRR>_2.fastq.gz` for paired-end libraries. The table in `fetchngs/output/samplesheet/samplesheet.csv` tells which GEO sample each file belongs to.

## Running crepas

```bash
mkdir -p ${HOME}/projects/scarseq_tutorial/crepas/
```

### Input samplesheet CSV

The samplesheet contains one row per library. Each SCAR-seq sample (`antibody` set) points to its stranded input with the `input_control` and `input_control_biological_replicate` columns, while the input rows leave `antibody` and `input_control` empty. Some of the libraries are single-end and others paired-end, so `fastq_2` is only filled for the latter. We save it as `scarseq_tutorial.crepas.samplesheet.csv`:

<details>
<summary>Show the samplesheet (56 rows)</summary>

```csv
sample,biological_replicate,input_control,input_control_biological_replicate,exp_type,strandedness,antibody,fastq_1,fastq_2
SCAR_410_T0_H3K27me3,1,SCAR_410_T0_H3K27me3_strandedInput,1,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725306_SRR12214703.fastq.gz,
SCAR_410_T0_H3K27me3,2,SCAR_410_T0_H3K27me3_strandedInput,2,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725315_SRR12214712.fastq.gz,
SCAR_410_T0_H3K27me3_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725342_SRR12214739.fastq.gz,
SCAR_410_T0_H3K27me3_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034080_SRR17874677_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034080_SRR17874677_2.fastq.gz
SCAR_410_T3_H3K27me3,1,SCAR_410_T3_H3K27me3_strandedInput,1,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725307_SRR12214704.fastq.gz,
SCAR_410_T3_H3K27me3,2,SCAR_410_T3_H3K27me3_strandedInput,2,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725316_SRR12214713.fastq.gz,
SCAR_410_T3_H3K27me3_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725343_SRR12214740.fastq.gz,
SCAR_410_T3_H3K27me3_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034081_SRR17874678_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034081_SRR17874678_2.fastq.gz
SCAR_410_T8_H3K27me3,1,SCAR_410_T8_H3K27me3_strandedInput,1,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725308_SRR12214705.fastq.gz,
SCAR_410_T8_H3K27me3,2,SCAR_410_T8_H3K27me3_strandedInput,2,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725317_SRR12214714.fastq.gz,
SCAR_410_T8_H3K27me3_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725344_SRR12214741.fastq.gz,
SCAR_410_T8_H3K27me3_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034082_SRR17874679_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034082_SRR17874679_2.fastq.gz
SCAR_439_T0_H3K27me3,1,SCAR_439_T0_H3K27me3_strandedInput,1,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725312_SRR12214709.fastq.gz,
SCAR_439_T0_H3K27me3,2,SCAR_439_T0_H3K27me3_strandedInput,2,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725321_SRR12214718.fastq.gz,
SCAR_439_T0_H3K27me3_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725348_SRR12214745.fastq.gz,
SCAR_439_T0_H3K27me3_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034086_SRR17874683_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034086_SRR17874683_2.fastq.gz
SCAR_439_T3_H3K27me3,1,SCAR_439_T3_H3K27me3_strandedInput,1,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725313_SRR12214710.fastq.gz,
SCAR_439_T3_H3K27me3,2,SCAR_439_T3_H3K27me3_strandedInput,2,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725322_SRR12214719.fastq.gz,
SCAR_439_T3_H3K27me3_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725349_SRR12214746.fastq.gz,
SCAR_439_T3_H3K27me3_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034087_SRR17874684_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034087_SRR17874684_2.fastq.gz
SCAR_439_T8_H3K27me3,1,SCAR_439_T8_H3K27me3_strandedInput,1,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725314_SRR12214711.fastq.gz,
SCAR_439_T8_H3K27me3,2,SCAR_439_T8_H3K27me3_strandedInput,2,SCAR-seq,forward,H3K27me3,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725323_SRR12214720.fastq.gz,
SCAR_439_T8_H3K27me3_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX8725350_SRR12214747.fastq.gz,
SCAR_439_T8_H3K27me3_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034088_SRR17874685_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034088_SRR17874685_2.fastq.gz
SCAR_410_T0_H3K27ac,1,SCAR_410_T0_H3K27ac_strandedInput,1,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14033999_SRR17874620_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14033999_SRR17874620_2.fastq.gz
SCAR_410_T0_H3K27ac,2,SCAR_410_T0_H3K27ac_strandedInput,2,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034059_SRR17874632_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034059_SRR17874632_2.fastq.gz
SCAR_410_T0_H3K27ac_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034071_SRR17874644_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034071_SRR17874644_2.fastq.gz
SCAR_410_T0_H3K27ac_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034011_SRR17874656_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034011_SRR17874656_2.fastq.gz
SCAR_410_T1_H3K27ac,1,SCAR_410_T1_H3K27ac_strandedInput,1,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034000_SRR17874621_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034000_SRR17874621_2.fastq.gz
SCAR_410_T1_H3K27ac,2,SCAR_410_T1_H3K27ac_strandedInput,2,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034060_SRR17874633_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034060_SRR17874633_2.fastq.gz
SCAR_410_T1_H3K27ac_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034072_SRR17874645_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034072_SRR17874645_2.fastq.gz
SCAR_410_T1_H3K27ac_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034012_SRR17874657_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034012_SRR17874657_2.fastq.gz
SCAR_410_T3_H3K27ac,1,SCAR_410_T3_H3K27ac_strandedInput,1,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034001_SRR17874622_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034001_SRR17874622_2.fastq.gz
SCAR_410_T3_H3K27ac,2,SCAR_410_T3_H3K27ac_strandedInput,2,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034061_SRR17874634_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034061_SRR17874634_2.fastq.gz
SCAR_410_T3_H3K27ac_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034073_SRR17874646_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034073_SRR17874646_2.fastq.gz
SCAR_410_T3_H3K27ac_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034013_SRR17874658_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034013_SRR17874658_2.fastq.gz
SCAR_410_T8_H3K27ac,1,SCAR_410_T8_H3K27ac_strandedInput,1,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034050_SRR17874623_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034050_SRR17874623_2.fastq.gz
SCAR_410_T8_H3K27ac,2,SCAR_410_T8_H3K27ac_strandedInput,2,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034062_SRR17874635_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034062_SRR17874635_2.fastq.gz
SCAR_410_T8_H3K27ac_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034002_SRR17874647_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034002_SRR17874647_2.fastq.gz
SCAR_410_T8_H3K27ac_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034014_SRR17874659_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034014_SRR17874659_2.fastq.gz
SCAR_439_T0_H3K27ac,1,SCAR_439_T0_H3K27ac_strandedInput,1,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034051_SRR17874624_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034051_SRR17874624_2.fastq.gz
SCAR_439_T0_H3K27ac,2,SCAR_439_T0_H3K27ac_strandedInput,2,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034063_SRR17874636_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034063_SRR17874636_2.fastq.gz
SCAR_439_T0_H3K27ac_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034003_SRR17874648_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034003_SRR17874648_2.fastq.gz
SCAR_439_T0_H3K27ac_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034015_SRR17874660_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034015_SRR17874660_2.fastq.gz
SCAR_439_T1_H3K27ac,1,SCAR_439_T1_H3K27ac_strandedInput,1,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034052_SRR17874625_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034052_SRR17874625_2.fastq.gz
SCAR_439_T1_H3K27ac,2,SCAR_439_T1_H3K27ac_strandedInput,2,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034064_SRR17874637_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034064_SRR17874637_2.fastq.gz
SCAR_439_T1_H3K27ac_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034004_SRR17874649_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034004_SRR17874649_2.fastq.gz
SCAR_439_T1_H3K27ac_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034016_SRR17874661_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034016_SRR17874661_2.fastq.gz
SCAR_439_T3_H3K27ac,1,SCAR_439_T3_H3K27ac_strandedInput,1,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034053_SRR17874626_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034053_SRR17874626_2.fastq.gz
SCAR_439_T3_H3K27ac,2,SCAR_439_T3_H3K27ac_strandedInput,2,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034065_SRR17874638_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034065_SRR17874638_2.fastq.gz
SCAR_439_T3_H3K27ac_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034005_SRR17874650_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034005_SRR17874650_2.fastq.gz
SCAR_439_T3_H3K27ac_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034017_SRR17874662_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034017_SRR17874662_2.fastq.gz
SCAR_439_T8_H3K27ac,1,SCAR_439_T8_H3K27ac_strandedInput,1,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034054_SRR17874627_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034054_SRR17874627_2.fastq.gz
SCAR_439_T8_H3K27ac,2,SCAR_439_T8_H3K27ac_strandedInput,2,SCAR-seq,forward,H3K27ac,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034066_SRR17874639_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034066_SRR17874639_2.fastq.gz
SCAR_439_T8_H3K27ac_strandedInput,1,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034006_SRR17874651_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034006_SRR17874651_2.fastq.gz
SCAR_439_T8_H3K27ac_strandedInput,2,,,SCAR-seq,forward,,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034018_SRR17874663_1.fastq.gz,<PATH>/projects/scarseq_tutorial/fetchngs/output/fastq/SRX14034018_SRR17874663_2.fastq.gz
```

</details>

> [!NOTE]
> `410` and `439` in the sample names are the two cell lines of the article (wild type and MCM2-2A), and `T0`, `T3` and `T8` the hours after the nascent DNA labelling, as in the GEO sample names.

### Parameter YAML file

The `okseq_rfd_file` and `initiation_zones` parameters point to the outputs of the [OK-seq tutorial](okseq.md). We also enable the counting of transposable elements (TEs), which needs a TE annotation in GTF format (`te_gtf`):

```yaml
input: <PATH>/projects/scarseq_tutorial/crepas/scarseq_tutorial.crepas.samplesheet.csv
outdir: <PATH>/projects/scarseq_tutorial/crepas/output
multiqc_title: SCAR-seq_H3K27me3_H3K27ac
genome: GRCm39
aligner: chromap
fasta: <GRCm39_FASTA>
gtf: <GRCm39_GTF>
blacklist: <GRCm39_BLACKLIST_BED>
te_gtf: <GRCm39_TE_GTF>
okseq_rfd_file: <PATH>/projects/okseq_tutorial/crepas/output/chromap/mergedLibrary/OK-seq/6_collect/OK-seq_SRR7535256_bRep_1.mLb.mkD.flT1.flT3.flTbl.rpm.RFD.w1000.s30.d30.z1.collect.flT_by_counts.tsv
initiation_zones: <PATH>/projects/okseq_tutorial/crepas/output/chromap/mergedLibrary/OK-seq/7_rfd_to_iz/OK-seq_SRR7535256_bRep_1.mLb.mkD.flT1.flT3.flTbl.rpm.RFD.w1000.s30.d30.z1.collect.flT_by_counts.init_zones.rm_overlaps.bed
seq_platform: ILLUMINA
read_length: 100
fragment_size: 1000
coverage_extend_reads: 1000
coverage_bin_size: 1000
map_n_multimappers: 100
save_align_intermeds: true
save_partition_intermeds: true
skip_genes_plotprofile: true
skip_plot_fingerprint: true
skip_te_counting: false
skip_telocal: false
```

> [!TIP]
> Check the exact names of the OK-seq files in your `okseq_tutorial/crepas/output/` directory: they include the filters and smoothing parameters used in that run.

Save this file as `scarseq_tutorial.crepas.params.yaml` in the `crepas/` subdirectory.

### SBATCH script

The SBATCH script is the one of the [OK-seq tutorial](okseq.md), replacing `okseq_tutorial` with `scarseq_tutorial`:

```bash
cd ${HOME}/projects/scarseq_tutorial/crepas/

nextflow run grothlab/crepas \
    -r 1.2.0 \
    -profile singularity \
    -params-file ${HOME}/projects/scarseq_tutorial/crepas/scarseq_tutorial.crepas.params.yaml
```

Save it as `scarseq_tutorial.crepas.sbatch.sh`.

### Submitting the job

```bash
cd ${HOME}/projects/scarseq_tutorial/crepas/
sbatch scarseq_tutorial.crepas.sbatch.sh
```

The partition tracks and tables are in `crepas/output/chromap/mergedLibrary/SCAR-seq/6_collect/`, the averages across biological replicates in its `bRep_average/` subdirectory, and the partition plots around the OK-seq initiation zones in `crepas/output/chromap/mergedLibrary/SCAR-seq/7_plots/`.
