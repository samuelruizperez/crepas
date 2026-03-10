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
    tuple val(meta), path("*.bed"), emit: bed
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.filtered"
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    def idr_threshold = args2.contains('--min_signal ') ? args2.split('--min_signal ')[1].split(' ')[0].toBigDecimal() : 0
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    IDR_THRESH_TRANSFORMED=\$(awk -v p=${idr_threshold} 'BEGIN{print -log(p)/log(10)}')

    awk \\
        ${args} \\
        'BEGIN{OFS="\\t"} \$12>=\${IDR_THRESH_TRANSFORMED} {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' \\
        ${peak_file} \\
        | sort \\
        | uniq \\
        | sort -k7n,7n \\
        > ${prefix}.${peak_type}


        '{if (\$3 != -1) print \$0}' \\
        $bed \\
        > ${prefix}.bed

    
IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')

awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${IDR_OUTPUT} | sort | uniq | sort -k7n,7n | gzip -nc > ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
