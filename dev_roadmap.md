
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



# FileAlreadyExistsException when using copyTo on a directory that already exists #3887

Sometimes, running the pipeline several times with -resume will cause the pipeline to fail with a FileAlreadyExistsException when using copyTo on a directory that already exists. This is because the copyTo method does not overwrite existing directories. See: https://github.com/nextflow-io/nextflow/discussions/3887#discussioncomment-5667052


# make BAM_CREATE_SCAR_PARTITIONS more readable