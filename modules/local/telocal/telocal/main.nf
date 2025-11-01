process TELOCAL {
    tag "$meta.id"
    label 'process_high_memory_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mhammelllab/telocal:latest' :
        'docker://mhammelllab/telocal:latest' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(genic_gtf_or_index)
    tuple val(meta3), path(te_gtf_or_index)
    val skip_gz

    output:
    tuple val(meta), path("*.cntTable*"),                        emit: counts
    path "versions.yml",                                        emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TElocal: \$(TElocal --version | sed 's/TElocal //g')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip[[:space:]]*//' )

    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd = skip_gz ? "touch ${prefix}.cntTable" : "touch ${prefix}.cntTable.gz"

    """
    ${cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TElocal: \$(TElocal --version | sed 's/TElocal //g')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip[[:space:]]*//' )
    END_VERSIONS
    """
}
