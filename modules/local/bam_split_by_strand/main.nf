/*
 * Split a BAM file by strand
 */
process BAM_SPLIT_BY_STRAND {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.forward.bam"), emit: f_bam
    tuple val(meta), path("*.reverse.bam"), emit: r_bam
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix            = task.ext.prefix ?: "${meta.id}"
    def strand_extension1 = meta.strandedness == 'reverse' ? 'reverse' : 'forward'
    def strand_extension2 = meta.strandedness == 'reverse' ? 'forward' : 'reverse'

    """
    samtools view \\
        -F 20 -h $bam | \\
            samtools view -Sb -h \\
                > ${prefix}.${strand_extension1}.bam

    samtools view \\
        -f 16 -h $bam | \\
            samtools view -Sb -h \\
                > ${prefix}.${strand_extension2}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
