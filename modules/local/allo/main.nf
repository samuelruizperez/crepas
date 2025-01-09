process ALLO {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/allo:1.2.0--pyhdfd78af_0' :
        'quay.io/biocontainers/allo:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam")             , emit: bam
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq = meta.single_end ? '-seq "se"' : '-seq "pe"'
    def random = meta.control ? '--random' : ''
    """
    allo \\
        $args \\
        $bam \\
        $seq \\
        $random \\
        -p $task.cpus \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        allo: \$(echo \$(allo -v 2>&1))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        allo: \$(echo \$(allo -v 2>&1))
    END_VERSIONS
    """
}
