process NORMALIZE_STRANDS {
    tag "$archive"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(tab), val(cpm)

    output:
    tuple val(meta), path("*.tab"), emit: tab
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk \\
        $args \\
        -v c=\${$cpm} \\
        'BEGIN{OFS="\t"}{print \$1, \$2, \$3, (\$4+1)/c, \$5, \$6}' \\
        $tab \\
        > ${prefix}.cpm.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.cpm.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | sed -e "s/awk v//g")
    END_VERSIONS
    """
}
