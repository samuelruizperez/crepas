/*
 * Split a BAM file by strand
 * It uses SAM FLAG 16: read reverse strand (0x10)
 */
process BAM_SPLIT_BY_STRAND {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.forward.bam"), emit: f_bam
    tuple val(meta), path("*.reverse.bam"), emit: r_bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix            = task.ext.prefix ?: "${meta.id}"
    
    if (meta.strandedness == 'forward') {
    """
        samtools view \\
            --threads ${task.cpus-1} \\
            --exclude-flags 20 \\
            --with-header \\
            --bam \\
            --output ${prefix}.forward.bam \\
            ${bam}

        samtools view \\
            --threads ${task.cpus-1} \\
            --require-flags 16 \\
            --with-header \\
            --bam \\
            --output ${prefix}.reverse.bam \\
            ${bam}
        """
    } else if (meta.strandedness == 'reverse') {
    """
        samtools view \\
            --threads ${task.cpus-1} \\
            --exclude-flags 20 \\
            --with-header \\
            --bam \\
            --output ${prefix}.reverse.bam \\
            ${bam}

        samtools view \\
            --threads ${task.cpus-1} \\
            --require-flags 16 \\
            --with-header \\
            --bam \\
            --output ${prefix}.forward.bam \\
            ${bam}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.forward.bam
    touch ${prefix}.reverse.bam
    """
}
