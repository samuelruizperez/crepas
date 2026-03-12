process FILTER_IDR_PEAKS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(peak_file)
    val peak_type
    val idr_threshold

    output:
    tuple val(meta), path("${prefix}.${peak_type}"), emit: filtered_peaks
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.filtered_idr"
    def idr_threshold_arg = idr_threshold ?: 0.05
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    IDR_THRESH_TRANSFORMED=\$(awk -v p=${idr_threshold_arg} 'BEGIN{print -log(p)/log(10)}')

    awk \\
        ${args} \\
        'BEGIN{OFS="\\t"} \$12>=\${IDR_THRESH_TRANSFORMED} {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' \\
        ${peak_file} \\
        | sort \\
        | uniq \\
        | sort -k7n,7n \\
        > ${prefix}.${peak_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.filtered_idr"
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    touch  ${prefix}.${peak_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
