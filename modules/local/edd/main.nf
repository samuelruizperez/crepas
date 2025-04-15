
process EDD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/edd:1.1.19--py27h4329609_4':
        'biocontainers/edd:1.1.19--py27h4329609_4' }"

    input:
    tuple val(meta), path(ip_bam), path(input_bam)
    tuple val(meta2), path(chrom_sizes)
    tuple val(meta3), path(unalignable_regions)

    output:
    tuple val(meta), path("*_peaks.bed"),                               emit: peaks
    tuple val(meta), path("*.log"),                                     emit: log
    tuple val(meta), path("*_bin_score.bedgraph"),      optional:true,  emit: bin_score
    tuple val(meta), path("*_log_ratio*.bedgraph"),     optional:true,  emit: log_ratios
    path  "versions.yml",                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    edd \\
        $args \\
        --nprocs $task.cpus \\
        --ip_bam $ip_bam \\
        --input_bam $input_bam \\
        --chrom_sizes $chrom_sizes \\
        --unalignable_regions $unalignable_regions \\
        --output_dir ./${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        edd: \$(edd --version | sed -e "s/edd //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        edd: \$(edd --version | sed -e "s/edd //g")
    END_VERSIONS
    """
}
