process ALLO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/allo_samtools:9c6e229f802e6c51' :
        'community.wave.seqera.io/library/allo_samtools:0cfc6883c80e7505' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam")      , emit: bam
    tuple val(meta), path("*.sam")      , emit: sam
    tuple val(meta), path("*.log")      , emit: log
    tuple val("${task.process}"), val('allo'), eval("allo -v 2>&1"), emit: versions_allo, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def seq     = meta.single_end ? '-seq "se"' : '-seq "pe"'
    """
    allo \\
        ${args} \\
        ${bam} \\
        ${seq} \\
        -p ${task.cpus} \\
        -o ${prefix}.sam \\
        &> >(tee ${prefix}.log >&2)

    samtools view \\
        -bS ${prefix}.sam \\
        > ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.sam
    touch  ${prefix}.bam
    touch  ${prefix}.log
    """
}
