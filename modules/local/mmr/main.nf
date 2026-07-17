process MMR {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/mmr:c5ce80a--0042edeccacdf2bc' :
        'community.wave.seqera.io/library/mmr:c5ce80a--1ec1f5037b507cc0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam")      , emit: bam
    tuple val("${task.process}"), val('mmr'), eval("echo $VERSION"), topic: versions, emit: versions_mmr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def pair_usage_arg     = meta.single_end ? "" : "--pair-usage"
    VERSION = 'c5ce80a' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mmr \\
        $args \\
        $pair_usage_arg \\
        --threads $task.cpus \\
        --verbose \\
        -o ${prefix}.bam \\
        $bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = 'c5ce80a' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch  ${prefix}.bam
    """
}
