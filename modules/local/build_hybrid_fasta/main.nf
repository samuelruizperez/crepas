process BUILD_HYBRID_FASTA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:22.04'
        : 'nf-core/ubuntu:22.04'}"

    input:
    tuple val(meta), path(fasta, stageAs: "genome/*"), val(genome)
    tuple val(meta2), path(spikein_fasta, stageAs: "spikein_genome/*"), val(spikein_genome)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    {
        cat "${fasta}"
        awk -v suffix="_${spikein_genome}" '/^>/{print \$0 suffix; next} 1' "${spikein_fasta}"
    } > "${genome}_${spikein_genome}.fa"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${genome}_${spikein_genome}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
