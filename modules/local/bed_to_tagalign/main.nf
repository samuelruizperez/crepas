//
// Converting BED to TAGALIGN format
// Based on: https://github.com/ENCODE-DCC/chip-seq-pipeline2; step 2a
//
process BED_TO_TAGALIGN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.tagAlign"), emit: tagalign
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
        awk \\
            ${args} \\
            'BEGIN{OFS="\\t"}{\$4="N";\$5="1000";print \$0}' \\
            ${bed} \\
            > ${prefix}.tagAlign

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
        END_VERSIONS
        """
    } else {
        """
        awk \\
            ${args} \\
            'BEGIN{OFS="\\t"}{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' \\
            ${bed} \\
            > ${prefix}.tagAlign

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.tagAlign

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
