
1. Handling of scaffolds for exclusion in partition plot generation

2. ratio SCAR-seq input subtraction

3. samtools stat summary table

5. peak analysis for SCAR-seq

samtools collate after resume generates error in allo cause it inputs all the temp .bams

Note:

igv and multiqc directories are created per allocation_method and per peak type

TODO: for ChOR-seq, what happens in MACS3 when the inputs have very variable number of reads due to issues with ligation/library construction.

ADD HISAT2 to the pipeline


- add consensus peak calling for ATAC-seq (see https://github.com/grothlab/glseq/blob/c753120ee33a5e0b5e8bb3dea319aa45eed34473/workflows/glseq.nf#L589)

and https://github.com/nf-core/atacseq/blob/dev/subworkflows/local/bed_consensus_quantify_qc_bedtools_featurecounts_deseq2.nf

- add peak calling when there are no controls


# FileAlreadyExistsException when using copyTo on a directory that already exists #3887

Sometimes, running the pipeline several times with -resume will cause the pipeline to fail with a FileAlreadyExistsException when using copyTo on a directory that already exists. This is because the copyTo method does not overwrite existing directories. See: https://github.com/nextflow-io/nextflow/discussions/3887#discussioncomment-5667052


# make BAM_CREATE_SCAR_PARTITIONS more readable

# update nf-validation to nf-schema


# for chorseq peak calling, check that https://www.nature.com/articles/s41596-021-00585-3#Sec60 is followed correctly, step 163 says to use the BAMs from step 162, which are not BAMs, but bedgraphs

# make sure only unique reads are kept in bams for chorseq rrpm calculation


# max_time etc. paramteres in output -params-file does not match schema pattern

--max_time: string [{days=3, hours=72, seconds=259200, millis=259200000, durationInMillis=259200000, minutes=4320}] does not match pattern ^(\d+\.?\s*(s|m|h|d|day)\s*)+$ ({days=3, hours=72, seconds=259200, millis=259200000, durationInMillis=259200000, minutes=4320})
* --max_memory: string [{mega=65536, kilo=67108864, giga=64, bytes=68719476736}] does not match pattern ^\d+(\.\d+)?\.?\s*(K|M|G|T)?B$ ({mega=65536, kilo=67108864, giga=64, bytes=68719476736})

# downsampling by min in type (chip or input) 
should be done
separately for each experiment (chipseq, scarseq, etc.) and not for all experiments together

# remove redundant samtools stats before downsampling

# add https://github.com/LHentges/LanceOtron

# add flt3 of spike in bams before spike-in normalization

# what happens with process_stats_summary.R if there is a dot in either the genome or the spikein_genome name?

# verify that befor SOI both files have same number of bins in the same order

# add input_basedir
See https://nf-co.re/pixelator/latest/docs/usage/