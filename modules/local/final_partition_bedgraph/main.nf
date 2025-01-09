process FINAL_PARTITION_BEDGRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.tmp"), emit: tmp
    tuple val(meta), path("*.txt"), emit: txt
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk \\
        $args \\
        '\$8>0 || \$14>0' \\
        $file \\
        | cut -f -3,8,14,20,26,30- \\
        > ${prefix}.txt

    awk \\
        $args2 \\
        '{printf "%s\\t%d\\t%d\\t%2.3f\\n" , \$1,\$2,\$3,\$9}' \\
        ${prefix}.txt \\
        > ${prefix}.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.txt
    touch  ${prefix}.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
