process IDR_FILTER_THRESHOLD {
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
    tuple val(meta), path("${prefix}.${peak_type}"), emit: peaks
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

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
        | LC_COLLATE=C sort -T '.' \\
        | uniq \\
        | sort -k7n,7n \\
        > ${prefix}.${peak_type}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.filtered_idr"
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    touch  ${prefix}.${peak_type}
    """
}
