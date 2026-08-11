/*
 * Split a BAM file by strand
 */
process BAM_SPLIT_BY_STRAND {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), val(exp_type), val(strandedness)

    output:
    tuple val(meta), path("*.forward.bam"), emit: f_bam
    tuple val(meta), path("*.reverse.bam"), emit: r_bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix            = task.ext.prefix ?: "${meta.id}"

    if (exp_type in ['SCAR-seq', 'OK-seq'] && strandedness == 'forward') {
    """
        # Forward strand: reads aligned in the forward orientation
        # Equivalent to: !4 && !16
        samtools view \\
            --threads ${task.cpus-1} \\
            --expr '!flag.unmap && !flag.reverse' \\
            --with-header \\
            --bam \\
            --output ${prefix}.forward.bam \\
            ${bam}

        # Reverse strand: reads aligned in the reverse orientation
        # Equivalent to: !4 && 16
        samtools view \\
            --threads ${task.cpus-1} \\
            --expr '!flag.unmap && flag.reverse' \\
            --with-header \\
            --bam \\
            --output ${prefix}.reverse.bam \\
            ${bam}
        """
    } else if (exp_type in ['SCAR-seq', 'OK-seq'] && strandedness == 'reverse') {
    """
        # Reverse strand: the insert strandedness is flipped, so reads aligned in the forward
        # orientation correspond to the reverse strand
        # Equivalent to: !4 && !16
        samtools view \\
            --threads ${task.cpus-1} \\
            --expr '!flag.unmap && !flag.reverse' \\
            --with-header \\
            --bam \\
            --output ${prefix}.reverse.bam \\
            ${bam}

        # Forward strand: reads aligned in the reverse orientation
        # Equivalent to: !4 && 16
        samtools view \\
            --threads ${task.cpus-1} \\
            --expr '!flag.unmap && flag.reverse' \\
            --with-header \\
            --bam \\
            --output ${prefix}.forward.bam \\
            ${bam}
        """
    } else if (exp_type == 'eSPAN' && !meta.single_end) {
    """
        # Forward strand: both mates of every pair whose READ2 aligns in the forward orientation
        # Equivalent to: !4 && 1 && ((128 && !16) || (64 && !32))
        # NOTE: READ2 aligns forward, but in eSPAN this pair represents the genomic reverse/minus/Crick nascent strand.
        samtools view \\
            --threads ${task.cpus-1} \\
            --expr '!flag.unmap && flag.paired && ((flag.read2 && !flag.reverse) || (flag.read1 && !flag.mreverse))' \\
            --with-header \\
            --bam \\
            --output ${prefix}.reverse.bam \\
            ${bam}

        # Reverse strand: both mates of every pair whose READ2 aligns in the reverse orientation
        # Equivalent to: !4 && 1 && ((128 && 16) || (64 && 32))
        # NOTE: READ2 aligns reverse, but in eSPAN this pair represents the genomic forward/plus/Watson nascent strand.
        samtools view \\
            --threads ${task.cpus-1} \\
            --expr '!flag.unmap && flag.paired && ((flag.read2 && flag.reverse) || (flag.read1 && flag.mreverse))' \\
            --with-header \\
            --bam \\
            --output ${prefix}.forward.bam \\
            ${bam}
        """
    } else if (exp_type == 'eSPAN' && meta.single_end) {
    """
        # Equivalent to: !4 && !1 && 16
        samtools view \\
            --threads ${task.cpus-1} \\
            --expr '!flag.unmap && !flag.paired && flag.reverse' \\
            --with-header \\
            --bam \\
            --output ${prefix}.reverse.bam \\
            ${bam}

        # Equivalent to: !4 && !1 && !16
        samtools view \\
            --threads ${task.cpus-1} \\
            --expr '!flag.unmap && !flag.paired && !flag.reverse' \\
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
