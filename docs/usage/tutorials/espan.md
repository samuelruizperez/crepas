# eSPAN

For this example, we will analyze the H3K4me3 eSPAN data of fission yeast (_Schizosaccharomyces pombe_) from the following article:

> Fang, Y., et al. (2024). Coordination of histone chaperones for parental histone segregation and epigenetic inheritance. _Genes & Development_, _38_(3-4), 189–204. [https://doi.org/10.1101/gad.351278.123](https://doi.org/10.1101/gad.351278.123)

The article compares wild-type cells with _mcm2-2A_, _dpb4Δ_ and _mrc1Δ_ mutants and their combinations. Each eSPAN library has a BrdU-IP-ssSeq library of the same cells as input control.

We start by creating a directory for this tutorial:

```bash
mkdir -p ${HOME}/projects/espan_tutorial
```

## Downloading the data

```bash
mkdir -p ${HOME}/projects/espan_tutorial/fetchngs/
```

### Input IDs file

We list the 24 SRA experiments (SRX) of the eSPAN and BrdU-IP-ssSeq libraries:

<details>
<summary>Show the 24 accessions</summary>

```csv
SRX21064165
SRX21064163
SRX21064166
SRX21064164
SRX21064186
SRX21064184
SRX21064190
SRX21064188
SRX21064194
SRX21064192
SRX21282706
SRX21282704
SRX22130984
SRX22130966
SRX22130988
SRX22130986
SRX22130960
SRX22130958
SRX22130964
SRX22130962
SRX22130992
SRX22130990
SRX22130996
SRX22130994
```

</details>

And we save it as `espan_tutorial.fetchngs.input.csv` in the `fetchngs` subdirectory.

### SBATCH script

For these runs, we download the files with SRA Tools (`--download_method sratools`) instead of the default FTP download from ENA:

```bash
nextflow run nf-core/fetchngs \
    -r 1.12.0 \
    -profile singularity \
    --input ${HOME}/projects/espan_tutorial/fetchngs/espan_tutorial.fetchngs.input.csv \
    --outdir ${HOME}/projects/espan_tutorial/fetchngs/output/ \
    --download_method sratools
```

Save it as `espan_tutorial.fetchngs.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/espan_tutorial/fetchngs/
sbatch espan_tutorial.fetchngs.sbatch.sh
```

## Running crepas

```bash
mkdir -p ${HOME}/projects/espan_tutorial/crepas/
```

### Input samplesheet CSV

All libraries are paired-end. We save the samplesheet as `espan_tutorial.crepas.samplesheet.csv`:

<details>
<summary>Show the samplesheet (24 rows)</summary>

```csv
sample,biological_replicate,input_control,input_control_biological_replicate,exp_type,strandedness,antibody,fastq_1,fastq_2
WT,1,WT_BrdU-IP-ssSeq,1,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064165_SRR25321978_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064165_SRR25321978_2.fastq.gz
WT_BrdU-IP-ssSeq,1,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064163_SRR25321980_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064163_SRR25321980_2.fastq.gz
WT,2,WT_BrdU-IP-ssSeq,2,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064166_SRR25321977_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064166_SRR25321977_2.fastq.gz
WT_BrdU-IP-ssSeq,2,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064164_SRR25321979_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064164_SRR25321979_2.fastq.gz
mrc1D,1,mrc1D_BrdU-IP-ssSeq,1,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064186_SRR25321957_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064186_SRR25321957_2.fastq.gz
mrc1D_BrdU-IP-ssSeq,1,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064184_SRR25321959_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064184_SRR25321959_2.fastq.gz
mrc1D,2,mrc1D_BrdU-IP-ssSeq,2,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064190_SRR25321953_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064190_SRR25321953_2.fastq.gz
mrc1D_BrdU-IP-ssSeq,2,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064188_SRR25321955_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064188_SRR25321955_2.fastq.gz
mcm2.2A_mrc1D,1,mcm2.2A_mrc1D_BrdU-IP-ssSeq,1,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064194_SRR25321949_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064194_SRR25321949_2.fastq.gz
mcm2.2A_mrc1D_BrdU-IP-ssSeq,1,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064192_SRR25321951_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21064192_SRR25321951_2.fastq.gz
mcm2.2A_mrc1D,2,mcm2.2A_mrc1D_BrdU-IP-ssSeq,2,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21282706_SRR25553667_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21282706_SRR25553667_2.fastq.gz
mcm2.2A_mrc1D_BrdU-IP-ssSeq,2,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21282704_SRR25553669_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX21282704_SRR25553669_2.fastq.gz
mcm2.2A,1,mcm2.2A_BrdU-IP-ssSeq,1,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130984_SRR26426179_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130984_SRR26426179_2.fastq.gz
mcm2.2A_BrdU-IP-ssSeq,1,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130966_SRR26426197_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130966_SRR26426197_2.fastq.gz
mcm2.2A,2,mcm2.2A_BrdU-IP-ssSeq,2,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130988_SRR26426175_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130988_SRR26426175_2.fastq.gz
mcm2.2A_BrdU-IP-ssSeq,2,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130986_SRR26426177_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130986_SRR26426177_2.fastq.gz
dpb4D,1,dpb4D_BrdU-IP-ssSeq,1,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130960_SRR26426203_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130960_SRR26426203_2.fastq.gz
dpb4D_BrdU-IP-ssSeq,1,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130958_SRR26426205_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130958_SRR26426205_2.fastq.gz
dpb4D,2,dpb4D_BrdU-IP-ssSeq,2,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130964_SRR26426199_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130964_SRR26426199_2.fastq.gz
dpb4D_BrdU-IP-ssSeq,2,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130962_SRR26426201_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130962_SRR26426201_2.fastq.gz
mcm2.2A_dpb4D,1,mcm2.2A_dpb4D_BrdU-IP-ssSeq,1,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130992_SRR26426171_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130992_SRR26426171_2.fastq.gz
mcm2.2A_dpb4D_BrdU-IP-ssSeq,1,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130990_SRR26426173_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130990_SRR26426173_2.fastq.gz
mcm2.2A_dpb4D,2,mcm2.2A_dpb4D_BrdU-IP-ssSeq,2,eSPAN,forward,H3K4me3,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130996_SRR26426167_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130996_SRR26426167_2.fastq.gz
mcm2.2A_dpb4D_BrdU-IP-ssSeq,2,,,eSPAN,forward,,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130994_SRR26426169_1.fastq.gz,<PATH>/projects/espan_tutorial/fetchngs/output/fastq/SRX22130994_SRR26426169_2.fastq.gz
```

</details>

### Parameter YAML file

For fission yeast, we give the PomBase genome sequence and its annotation in GFF3 format (`gff`, instead of `gtf`). The eSPAN partitions are plotted around known replication origins (`initiation_zones`). Since the yeast genome is small, we also use smaller windows and bins than for mouse:

```yaml
input: <PATH>/projects/espan_tutorial/crepas/espan_tutorial.crepas.samplesheet.csv
outdir: <PATH>/projects/espan_tutorial/crepas/output
multiqc_title: eSPAN_yeast
genome: ASM294v3
aligner: chromap
fasta: <SPOMBE_FASTA>
gff: <SPOMBE_GFF>
initiation_zones: <SPOMBE_REPLICATION_ORIGINS_BED>
seq_platform: ILLUMINA
read_length: 50
fragment_size: 250
partition_fragment_size: 250
partition_window_size: 100
partition_step_size: 100
partition_plot_range: 10
partition_iz_rm_overlap_range: 0
coverage_extend_reads: 250
coverage_bin_size: 1
save_align_intermeds: true
save_partition_intermeds: true
skip_genes_plotprofile: true
skip_plot_fingerprint: true
skip_flT3: false
skip_flTbl: true
```

Save this file as `espan_tutorial.crepas.params.yaml` in the `crepas/` subdirectory.

### SBATCH script

```bash
cd ${HOME}/projects/espan_tutorial/crepas/

nextflow run grothlab/crepas \
    -r 1.2.0 \
    -profile singularity \
    -params-file ${HOME}/projects/espan_tutorial/crepas/espan_tutorial.crepas.params.yaml
```

Save it as `espan_tutorial.crepas.sbatch.sh`, with the same `#SBATCH` lines and modules as in the [OK-seq tutorial](okseq.md).

### Submitting the job

```bash
cd ${HOME}/projects/espan_tutorial/crepas/
sbatch espan_tutorial.crepas.sbatch.sh
```

The partition tracks and plots are in `crepas/output/chromap/mergedLibrary/eSPAN/`, as for SCAR-seq.
