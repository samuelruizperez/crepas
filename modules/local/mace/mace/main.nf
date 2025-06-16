
process MACE_MACE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mace:1.2--py27he7e273a_2':
        'biocontainers/mace:1.2_cv1' }"

    input:
    tuple val(meta), path(forward_bw), path(reverse_bw)
    tuple val(meta2), path(chrom_sizes)

    output:
    tuple val(meta), path("*.border.bed"),             emit: border
    tuple val(meta), path("*.border_cluster.bed"),     emit: border_cluster
    tuple val(meta), path("*.border_pair_elite.bed"),  emit: border_pair_elite
    tuple val(meta), path("*.border_pair.bed"),        emit: border_pair
    path  "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mace.py \\
        ${args} \\
        --chromSize ${chrom_sizes} \\
        --forward ${forward_bw} \\
        --reverse ${reverse_bw} \\
        --out-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mace: \$(mace.py --version | sed -e "s/mace //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.border.bed
    touch ${prefix}.border_cluster.bed
    touch ${prefix}.border_pair_elite.bed
    touch ${prefix}.border_pair.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mace: \$(mace.py --version | sed -e "s/mace //g")
    END_VERSIONS
    """
}
