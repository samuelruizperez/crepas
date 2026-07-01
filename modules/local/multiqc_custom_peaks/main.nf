process MULTIQC_CUSTOM_PEAKS {
    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/48492b71efcc81e0c02b7319ed42b67e4cf12db4350838e869157d0305a7c890/data' :
        'community.wave.seqera.io/library/coreutils_gawk:c6711ba4a1d2d075' }"

    input:
    tuple val(meta), path(peak), path(frip)
    path peak_count_header
    path frip_score_header

    output:
    tuple val(meta), path("*.peak_count_mqc.tsv"), emit: count
    tuple val(meta), path("*.FRiP_mqc.tsv")      , emit: frip
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $peak | wc -l | awk -v OFS='\t' '{ print "${prefix}", \$1 }' | cat $peak_count_header - > ${prefix}.peak_count_mqc.tsv
    cat $frip_score_header $frip > ${prefix}.FRiP_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
