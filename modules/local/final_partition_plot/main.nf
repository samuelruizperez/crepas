process FINAL_PARTITION_PLOT {
    tag "$meta.id"
    label 'process_medium_memory'

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda (params.enable_conda ? "conda-forge::r-base bioconda::bioconductor-deseq2 bioconda::bioconductor-biocparallel bioconda::bioconductor-tximport bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_r-argparse_r-dplyr_pruned:aad6cf5716386302' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_r-argparse_r-dplyr_pruned:ff63fd989740e4c5' }"

    input:
    tuple val(meta), path(partition), path(strandedinput), path(scarminusinput), path(okazaki)
    path blacklist
    path initiation_zones
    val scaffolds

    output:
    path "*_scatter_plots.pdf"                , optional:true, emit: scatter_pdf
    path "*_partition_RAW.pdf"                , emit: partition_raw_pdf
    path "*_partition_SMOOTHED.pdf"           , emit: partition_smoothed_pdf
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def okazaki_arg   = okazaki.isEmpty() ? '' : "--okazaki_file $okazaki"
    // if --exclude_chromosomes is in args and remove_scaffolds is true, then append scaffolds to the comma-separated list of --exclude_chromosomes
    def exclude_chromosomes = args.contains("--exclude_chromosomes") && scaffolds != null ? "--exclude_chromosomes ${args.split("--exclude_chromosomes")[1].split(" ")[0]},${scaffolds}" : args
    """
    SCAR_partition_plots.R \\
        --scar_partition_file $partition \\
        --scarminusinput_partition_file $scarminusinput \\
        --strandedinput_partition_file $strandedinput \\
        $okazaki_arg \\
        --initiation_zones $initiation_zones \\
        --blacklist $blacklist \\
        --prefix $prefix \\
        --outdir ./ \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
