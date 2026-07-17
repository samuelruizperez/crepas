process CHROMSIZES_SPLIT_BY_GENOME {
    tag "$sizes"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data' :
        'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9' }"

    input:
    tuple val(meta), path(sizes)
    val exo_genome_string
    val endo_genome_string

    output:
    tuple val(meta), path ("*.${endo_genome_string}.sizes") , emit: endo_sizes
    tuple val(meta), path ("*.${exo_genome_string}.sizes")  , emit: exo_sizes
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), topic: versions, emit: versions_awk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sizes.baseName}"
    """
    awk \\
        ${args} \\
        '!(\$1 ~ /_${exo_genome_string}\$/)' \\
        ${sizes} \\
        > ${prefix}.${endo_genome_string}.sizes

    awk \\
        ${args} \\
        '(\$1 ~ /_${exo_genome_string}\$/) {sub(/_${exo_genome_string}\$/, "", \$1); print}' \\
        ${sizes} \\
        > ${prefix}.${exo_genome_string}.sizes
    """

    stub:
    def prefix = task.ext.prefix ?: "${sizes.baseName}"
    """
    touch ${prefix}.${endo_genome_string}.sizes
    touch ${prefix}.${exo_genome_string}.sizes
    """
}
