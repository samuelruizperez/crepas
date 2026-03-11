process IDR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/idr:2.0.4.2--py39hcbe4a3b_5' :
        'biocontainers/idr:2.0.4.2--py39hcbe4a3b_5' }"

    input:
    tuple val(meta), path(sample_peaks), path(peak_list)
    val peak_type

    output:
    tuple val(meta), path("*.idrValues.txt"), emit: idr
    tuple val(meta), path("*.log.txt"      ), emit: log
    tuple val(meta), path("*.png"         ), emit: png
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.idr"
    def peak_list_arg = peak_list ? "--peak-list ${peak_list}" : ''
    if (sample_peaks.toList().size < 2) {
        log.error "[ERROR] idr needs at least two replicates only one provided."
    }
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    idr \\
        --samples ${sample_peaks} \\
        ${peak_list_arg} \\
        --input-file-type ${peak_type} \\
        --output-file ${prefix}.idrValues.txt \\
        --log-output-file ${prefix}.log.txt \\
        --plot \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        idr: \$(echo \$(idr --version 2>&1) | sed 's/^.*IDR //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.idr"
    if (sample_peaks.toList().size < 2) {
        log.error "[ERROR] idr needs at least two replicates only one provided."
    }
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    touch "${prefix}.idrValues.txt"
    touch "${prefix}.log.txt"
    touch "${prefix}.png"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        idr: \$(echo \$(idr --version 2>&1) | sed 's/^.*IDR //; s/ .*\$//')
    END_VERSIONS
    """
}
