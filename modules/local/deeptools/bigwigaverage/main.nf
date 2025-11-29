process DEEPTOOLS_BIGWIGAVERAGE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(bigwigs)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig, optional: true
    tuple val(meta), path("*.bedgraph"), emit: bedgraph, optional: true
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.bigwigAverage"
    def blacklist_cmd = blacklist ? "--blackListFileName ${blacklist}" : ""
    """
    bigwigAverage \\
        ${args} \\
        --bigwigs ${bigwigs.join(' ')} \\
        --numberOfProcessors ${task.cpus} \\
        --outFileName ${prefix}.bigWig \\
        $blacklist_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bigwigAverage --version | sed -e "s/bigwigAverage //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.bigwigAverage"
    """
    touch ${prefix}.bigwigAverage.bigWig
    touch ${prefix}.bigwigAverage.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bigwigAverage --version | sed -e "s/bigwigAverage //g")
    END_VERSIONS
    """
}