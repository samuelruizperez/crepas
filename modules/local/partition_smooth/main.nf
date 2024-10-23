process PARTITION_SMOOTH {
    tag "$archive"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(file)
    val   extension

    output:
    tuple val(meta), path("*.${extension}"), emit: sorted
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def buffer   = task.memory ? "--buffer-size=${task.memory.toGiga().intdiv(2)}G" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    LC_COLLATE=C sort \\
        $args \\
        $args2 \\
        --parallel=$task.cpus \\
        $buffer \\
        $file \\
        > ${prefix}.${extension}

    partition_smooth.pl \\
        ${TMP_DIR}/${CHR}_windows_F_CPM_minusinput.tab \\
        ${TMP_DIR}/${CHR}_windows_R_CPM_minusinput.tab \\
        ${RADIUS} ${DRADIUS} ${ZRADIUS} \\
        > ${TMP_DIR}/${CHR}_windows_RFD_minusinput.txt &

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(echo \$(sort --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed -e "s/sort v//g")
    END_VERSIONS
    """
}
