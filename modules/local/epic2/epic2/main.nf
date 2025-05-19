
process EPIC2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce0423a5369408ae7871e9988bb6d1e98539f1b1f5c6b60fa2f1d51d1a0bdcc5/data':
        'community.wave.seqera.io/library/epic2_pyranges:7b9f50ca2b9b5c13' }"

    input:
    tuple val(meta), path(treatment_bam), path(control_bam)
    tuple val(meta2), path(chrom_sizes)
    val effective_genome_fraction

    output:
    tuple val(meta), path("*.bed")                     , emit: bed
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chromsizes = chrom_sizes ? "--chromsizes $chrom_sizes" : ""
    def egf = effective_genome_fraction ? "--effective-genome-fraction $effective_genome_fraction" : ""
    def format = meta.single_end ? '' : args.contains('--guess-bampe') ? '' : '--guess-bampe'
    def treatment = treatment_bam   ? "--treatment ${treatment_bam.join(' ')}"   : ""
    def control = control_bam     ? "--control ${control_bam.join(' ')}" : ""
    """
    epic2 \\
        ${args} \\
        ${chromsizes} \\
        ${egf} \\
        ${treatment} \\
        ${control} \\
        --output ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epic2: \$(epic2 --version | sed -e "s/epic2 //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epic2: \$(epic2 --version | sed -e "s/epic2 //g")
    END_VERSIONS
    """
}
