process TELOCAL {
    tag "$meta.id"
    label 'process_high_memory_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mhammelllab/telocal:1.1.2' :
        'docker://mhammelllab/telocal:1.1.2' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(genic_gtf_or_index)
    tuple val(meta3), path(te_gtf_or_index)
    val skip_gz

    output:
    tuple val(meta), path("*.cntTable*"),                        emit: counts
    tuple val("${task.process}"), val('TElocal'), eval("TElocal --version 2>&1 | grep -Eo 'TElocal[[:space:]]+[0-9]+(\.[0-9]+)+' | tail -n 1 | sed 's/^TElocal[[:space:]]*//'"), emit: versions_telocal, topic: versions
    tuple val("${task.process}"), val('gzip'), eval("gzip --version | sed -n '1s/.*gzip[[:space:]]*//p'"), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args         ?: ""
    def prefix          = task.ext.prefix       ?: "${meta.id}"
    def gzip_cmd        = skip_gz ? "" : "gzip -f ${prefix}.cntTable"
    """
    TElocal \\
        ${args} \\
        --BAM ${bam} \\
        --GTF ${genic_gtf_or_index} \\
        --TE ${te_gtf_or_index} \\
        --project ${prefix}

    ${gzip_cmd}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd = skip_gz ? "touch ${prefix}.cntTable" : "touch ${prefix}.cntTable.gz"

    """
    ${cmd}
    """
}
