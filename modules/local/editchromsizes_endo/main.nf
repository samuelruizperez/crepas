process EDITCHROMSIZES_ENDO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(sizes)
    val filter_genome_string
    val keep_genome_string
    //val remove_scaffolds

    output:
    tuple val(meta), path ("*.sizes")   , emit: sizes
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    //def rm_command = remove_scaffolds ? " || \$1 ~ /\\./" : ''
    //'!(\$1 ~ /_${genome_string}\$/$rm_command)' \\
    """
    awk \\
        $args \\
        '!(\$1 ~ /_${filter_genome_string}\$/)' \\
        $sizes \\
        > ${sizes.baseName}.${keep_genome_string}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch endo.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
