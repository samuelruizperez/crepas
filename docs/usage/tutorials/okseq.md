# OK-seq

For this example, we will replicate the initiation zone and fork directionality analyses from the following article:

> Petryk, N., Dalby, M., Wenger, A., Stromme, C. B., Strandsby, A., Andersson, R., & Groth, A. (2018). MCM2 promotes symmetric inheritance of modified histones during DNA replication. _Science_, _361_(6409), 1389–1392. [https://doi.org/10.1126/science.aau0294](https://doi.org/10.1126/science.aau0294)

We start by creating a directory for this tutorial:

```bash
mkdir -p ${HOME}/projects/okseq_tutorial
```

## Downloading the data

First, we need to download the OK-seq data. We will do so using the [nf-core/fetchngs](https://nf-co.re/fetchngs/latest) pipeline, which works very well to fetch metadata and raw FastQ files from public databases, including the Gene Expression Omnibus (GEO).

We create a subdirectory for nf-core/fetchngs files:

```bash
mkdir -p ${HOME}/projects/okseq_tutorial/fetchngs/
```

### Input IDs file

We then create a CSV file with the SRA IDs of the samples we want to download.
In this case, we want the OK-seq data from the Petryk et al. 2018 paper, which is available under accession [SRR7535256](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR7535256&display=download). Therefore, the file CSV should look like this:

```csv
SRR7535256
```

And we save it as `okseq_tutorial.fetchngs.input.csv` in the `fetchngs` subdirectory.

### SBATCH script

Now, we will create an **SBATCH script** to submit the job to the cluster queue. You can use the following template.

```bash
#!/bin/bash

#SBATCH --job-name=FETCHNGS            # specify a name for the job
#SBATCH --mail-type=END,FAIL           # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=NONE               # email address to receive the notifications
#SBATCH -c 1                           # number of requested cores for the Nextflow head job
#SBATCH --mem=4gb                      # total requested RAM for the Nextflow head job
#SBATCH --time=2-00:00:00              # max. running time of the pipeline job, format in D-HH:MM:SS
#SBATCH --output=fetchngs_job.%j.log   # standard output and error log, '%j' gives the job ID

# Load the required modules
module load java/21.0.7 nextflow/25.10.4 singularity/3.8.7

# Create an output directory for the pipeline's output
mkdir -p ${HOME}/projects/okseq_tutorial/fetchngs/output/
cd ${HOME}/projects/okseq_tutorial/fetchngs/output/

# Run the pipeline
nextflow run nf-core/fetchngs \
    -r 1.12.0 \
    -profile singularity \
    --input ${HOME}/projects/okseq_tutorial/fetchngs/okseq_tutorial.fetchngs.input.csv \
    --outdir ${HOME}/projects/okseq_tutorial/fetchngs/output/
```

> [!TIP]
> On a cluster, you will usually also want Nextflow to submit each task as its own job instead of running everything inside the head job. Use your institution's [nf-core/configs](https://nf-co.re/configs) profile, if there is one, instead of (or in addition to) `singularity`.

Then place the script in your project directory, e.g., as `okseq_tutorial.fetchngs.sbatch.sh`.

Your folder structure should so far look something like this:

```
projects/
└── okseq_tutorial/
    └── fetchngs/
        ├── okseq_tutorial.fetchngs.input.csv
        └── okseq_tutorial.fetchngs.sbatch.sh
```

### Submitting the job

We can then submit the job to the queue with:

```bash
cd ${HOME}/projects/okseq_tutorial/fetchngs/
sbatch okseq_tutorial.fetchngs.sbatch.sh
```

You can then monitor the job status using `squeue` or by looking at the `fetchngs_job.%j.log` file in the `fetchngs` directory. You can also check the `.nextflow.log` file in the `fetchngs/output/` directory for more detailed logs.

When this job is finished, you should have the raw FastQ file `SRX4403092_SRR7535256.fastq.gz` downloaded in the `fetchngs/output/fastq/` subdirectory.

## Running crepas

Now that we have the data, we can run the _crepas_ pipeline to analyze it. We will create another subdirectory for the _crepas_ files:

```bash
mkdir -p ${HOME}/projects/okseq_tutorial/crepas/
```

### Input samplesheet CSV

Then, we need to create a samplesheet CSV file for the _crepas_ pipeline. This file should contain information about the samples we want to analyze. For this example, we will create a file called `okseq_tutorial.crepas.samplesheet.csv` in the `crepas` subdirectory with the following content:

```csv
sample,fastq_1,biological_replicate,exp_type,strandedness
SRR7535256,<PATH>/projects/okseq_tutorial/fetchngs/output/fastq/SRX4403092_SRR7535256.fastq.gz,1,OK-seq,forward
```

> [!NOTE]
> Make sure to replace `<PATH>` with the actual **full path** to your project directory.

### Parameter YAML file

Next, since the crepas pipeline has many more parameters, we will compile them into a parameter YAML file for easier management. The file should contain:

```yaml
input: <PATH>/projects/okseq_tutorial/crepas/okseq_tutorial.crepas.samplesheet.csv
outdir: <PATH>/projects/okseq_tutorial/crepas/output
multiqc_title: OK_seq_SRR7535256
genome: GRCm39
aligner: chromap
fasta: <GRCm39_FASTA>
gtf: <GRCm39_GTF>
blacklist: <GRCm39_BLACKLIST_BED>
seq_platform: ILLUMINA
read_length: 75
fragment_size: 200
partition_fragment_size: 200
save_partition_intermeds: true
trim_full_length: false
trim_minimum_length: 10
trim_maximum_length: 70
trim_nextseq: 0
trim_stringency: 1
trim_error_rate: 0.1
trim_q_cutoff: 20
skip_plot_fingerprint: true
skip_genes_plotprofile: true
skip_spp_qc: true
```

> [!NOTE]
> You can of course adjust any of these parameters or the other parameters of the pipeline, but in this case, since we want our results to be as similar as possible to those in the original Petryk et al. 2018 paper, we will set them as shown above.

Save this file as `okseq_tutorial.crepas.params.yaml` in the `crepas/` subdirectory.

### SBATCH script

Now, we will create an **SBATCH script** to submit the _crepas_ job to the cluster queue. You can use the following template.

```bash
#!/bin/bash

#SBATCH --job-name=CREPAS              # specify a name for the job
#SBATCH --mail-type=END,FAIL           # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=NONE               # email address to receive the notifications
#SBATCH -c 1                           # number of requested cores for the Nextflow head job
#SBATCH --mem=4gb                      # total requested RAM for the Nextflow head job
#SBATCH --time=4-00:00:00              # max. running time of the pipeline job, format in D-HH:MM:SS
#SBATCH --output=crepas_job.%j.log     # standard output and error log, '%j' gives the job ID

# Load the required modules
module load java/21.0.7 nextflow/25.10.4 singularity/3.8.7

# Run the pipeline from the crepas subdirectory, where Nextflow keeps its work/ directory and logs
cd ${HOME}/projects/okseq_tutorial/crepas/

nextflow run grothlab/crepas \
    -r 1.2.0 \
    -profile singularity \
    -params-file ${HOME}/projects/okseq_tutorial/crepas/okseq_tutorial.crepas.params.yaml
```

Save it as `okseq_tutorial.crepas.sbatch.sh` in the `crepas/` subdirectory. Your folder structure should now look like this:

```
projects/
└── okseq_tutorial/
    ├── fetchngs/
    │   ├── output/
    │   ├── okseq_tutorial.fetchngs.input.csv
    │   └── okseq_tutorial.fetchngs.sbatch.sh
    └── crepas/
        ├── okseq_tutorial.crepas.params.yaml
        ├── okseq_tutorial.crepas.samplesheet.csv
        └── okseq_tutorial.crepas.sbatch.sh
```

### Submitting the job

We can then submit the job to the queue with:

```bash
cd ${HOME}/projects/okseq_tutorial/crepas/
sbatch okseq_tutorial.crepas.sbatch.sh
```

As before, you can monitor the job with `squeue`, the `crepas_job.%j.log` file and the `.nextflow.log` file in the `crepas/` directory. If the run stops (for example, because it reached the time limit), add `-resume` to the `nextflow run` command and submit it again: Nextflow will reuse the tasks that already finished.

When the job is finished, the results are in `crepas/output/`, and the MultiQC report is in the `crepas/output/multiqc/` directory. The replication fork directionality (RFD) tracks are in `crepas/output/chromap/mergedLibrary/OK-seq/6_collect/`, and the initiation zones called from them in `crepas/output/chromap/mergedLibrary/OK-seq/7_rfd_to_iz/`. See the [output documentation](../../output.md) for a description of all output files.
