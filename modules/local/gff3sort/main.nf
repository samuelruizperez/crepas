process GFF3SORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7e888288455c539735f35c165c8b4f914aa986d283c239056bbdcbe932349875/data':
        'community.wave.seqera.io/library/gff3sort:0.1.a1a2bc9--3bdc1d11d91dce42' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}.sorted"
    """
    gff3sort.pl \\
        ${args} \\
        ${gtf} \\
        > ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}.sorted"
    """
    touch ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}
