/*
 * Split a BAM file by genome (in this case a string appended to chromosome names)
 */
process CHOR_NORM_FACTOR_CALCULATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_samtools_gawk:f298a1d70bb3367c' :
        'community.wave.seqera.io/library/htslib_samtools_gawk:f24fff5cf6f5ae5b' }"

    input:
    tuple val(meta), path(bam)
 
    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"
    """
    # Calculate normalization factor by dividing 10^6 by the total number of unique reads
    samtools view -c \\
        -F 260 $bam | \\
            awk '{print 1000000/\$1}' > ${prefix}.txt    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
