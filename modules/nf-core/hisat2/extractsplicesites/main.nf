process HISAT2_EXTRACTSPLICESITES {
    tag "$gtf"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c36472269e8898f63b7b65dd40433462d541f9e75f9401f0bf8488021275d006/data' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.splice_sites.txt"), emit: txt
    tuple val("${task.process}"), val('hisat2'), eval('hisat2 --version | grep -o "version [^ ]*" | cut -d " " -f 2'), topic: versions, emit: versions_hisat2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    hisat2_extract_splice_sites.py \\
        $args \\
        $gtf \\
        > ${gtf.baseName}.splice_sites.txt
    """

    stub:
    """
    touch ${gtf.baseName}.splice_sites.txt
    """
}
