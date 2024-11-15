process COLLECT_PARTITIONS_BY_CHROMOSOME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    paste \\
        $args \\
        ${files.join(' ')} \\
            > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paste: \$(echo \$(paste --version 2>&1) | sed 's/^.*(paste (GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paste: \$(echo \$(paste --version 2>&1) | sed 's/^.*(paste (GNU coreutils)) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
