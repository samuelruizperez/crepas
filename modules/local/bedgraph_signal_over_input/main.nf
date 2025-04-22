process BEDGRAPH_SIGNAL_OVER_INPUT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(ip_bedgraph), path(control_bedgraph)

    output:
    tuple val(meta), path("*.bedgraph") , emit: bedgraph
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    paste \\
        $args \\
        $ip_bedgraph \\
        $control_bedgraph \\
        | awk \\
            $args2 \\
            '{
                if (\$8 == 0) {
                    ratio = 0
                } else {
                    ratio = \$4 / \$8
                }
            print \$1, \$2, \$3, ratio
            }' OFS="\\t" \\
    > ${prefix}.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
