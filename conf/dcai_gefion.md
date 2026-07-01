# nf-core/configs: DCAI GEFION configuration

To use, run the pipeline with `-profile dcai_gefion`. This will download and launch the [`dcai_gefion.config`](../conf/dcai_gefion.config) which has been pre-configured with a setup suitable for the DCAI GEFION cluster.

## Running Nextflow workflow on DCAI GEFION

Various versions of Nextflow are available in the cluster as modules. You can see a list of available version using `module avail Nextflow` and then load your preferred Nextflow version using `module load`.

Nextflow shouldn't run directly on the login/submission node but on a compute node.

To do so make a shell script with a similar structure to the following code and submit with `sbatch my_script.sh`

```bash
#!/bin/bash

#SBATCH --job-name=<job_name>       # specify a name for the job
#SBATCH --mail-type=END,FAIL        # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=NONE            # email address to receive the notifications
#SBATCH -c 1                        # number of requested cores for the Nextflow head job (1 should be enough)
#SBATCH --mem=15gb                   # total requested RAM for the Nextflow head job (15 GB should be enough)
#SBATCH --time=0-02:00:00           # max. running time of the pipeline job, format in D-HH:MM:SS
#SBATCH --output=<job_name>.%j.log  # standard output and error log, '%j' gives the job ID
#SBATCH --account=<slurm_account>   # slurm account to submit this job
#SBATCH --reservation=<slurm_reservation>


# Set the SBATCH_ACCOUNT to your corresponding account in DCAI GEFION
export SBATCH_ACCOUNT=$(sacctmgr show association where users=$USER format=account -n -P)

# Set the SBATCH_RESERVATION if you have one
#export SBATCH_RESERVATION=<slurm_reservation>

# Set memory limits for the Nextflow head job
export NXF_OPTS="-Xms2g -Xmx4g"
export NXF_JVM_ARGS='-Xms2g -Xmx4g'

# Load the required modules
module purge
module load Apptainer/1.3.6 Nextflow/25.10.2

# Create an output directory for the pipeline run if it does not exist
mkdir -p <path_to_project_directory>/output/
cd <path_to_project_directory>/output/

# Run a public pipeline
nextflow run <pipeline_repo>/<pipeline_name> \
    -r <pipeline_version> \
    -profile dcai_gefion \
    -params-file <path_to_project_directory>/<params_file_yaml>
```