process SUBSTRACT_INPUT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(tab1), path(tab2), val(cpm)

    output:
    tuple val(updatedMeta), path("*.tab"), emit: tab
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    updatedMeta = meta + ['minusinput':true]  // Add minusinput to meta

    """
    paste \\
        $args \\
        $tab1 $tab2 \\
        | awk \\
            $args2 \\
            -v c=$cpm \\
            'BEGIN{OFS="\\t"}{if ((\$4-\$10) > 1/c) print \$1, \$2, \$3, \$4-\$10, \$5, \$6; else print \$1, \$2, \$3, 1/c, \$5, \$6}' - \\
            > ${prefix}.minusinput.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.minusinput.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
