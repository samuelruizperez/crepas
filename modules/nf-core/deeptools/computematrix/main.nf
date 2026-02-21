process DEEPTOOLS_COMPUTEMATRIX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bigwigs), path(beds)

    output:
    tuple val(meta), path("*.mat.gz") , emit: matrix
    tuple val(meta), path("*.mat.tab"), emit: table
    tuple val("${task.process}"), val('deeptools'), eval("computeMatrix --version | sed -e 's/computeMatrix //g'"), emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bigwigs_arg = "--scoreFileName ${bigwigs.join(' ')}"
    def beds_arg = "--regionsFileName ${beds.join(' ')}"

    """
    computeMatrix \\
        ${args} \\
        ${bigwigs_arg} \\
        ${beds_arg} \\
        --numberOfProcessors ${task.cpus} \\
        --outFileName ${prefix}.computeMatrix.mat.gz \\
        --outFileNameMatrix ${prefix}.computeMatrix.vals.mat.tab
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.computeMatrix.mat.gz
    touch ${prefix}.computeMatrix.vals.mat.tab
    """
}
