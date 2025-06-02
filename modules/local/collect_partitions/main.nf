//
// Collect windows, bigwigaverageoverbed scores, normalized scores, and RFD scores into a partition file
// Generate also a bedGraph file from the partition file
//
process COLLECT_PARTITIONS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(windows), path(bwaob_fwd), path(bwaob_rev), path(norm_or_smi_fwd), path(norm_or_smi_rev), path(rfd)

    output:
    tuple val(meta), path("*.tsv") , emit: tsv
    tuple val(meta), path("*.bdg") , emit: bdg
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def buffer   = task.memory ? "--buffer-size=${task.memory.toGiga().intdiv(2)}G" : ''
    def sort_cmd = "| LC_ALL=C sort --parallel=$task.cpus $buffer -k1,1 -k2,2n"

    """
    awk '
        FNR==NR { win[\$4]=\$1"\\t"\$2"\\t"\$3; next }
        FILENAME==ARGV[2] { bwaob_fwd[\$1]=\$4; next }
        FILENAME==ARGV[3] { bwaob_rev[\$1]=\$4; next }
        FILENAME==ARGV[4] { norm_fwd[\$1]=\$4; next }
        FILENAME==ARGV[5] { norm_rev[\$1]=\$4; next }
        FILENAME==ARGV[6] { rfd[\$1]=\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6; next }
        END {
            for (key in win) {
                if (bwaob_fwd[key] > 0 && bwaob_rev[key] > 0) {
                    print win[key] "\\t" bwaob_fwd[key] "\\t" bwaob_rev[key] "\\t" norm_fwd[key] "\\t" norm_rev[key] "\\t" rfd[key]
                }
            }
        }
    ' ${windows} ${bwaob_F} ${bwaob_R} ${norm_F} ${norm_R} ${rfd} \\
    ${sort_cmd} \\
    > ${prefix}.tsv

    # Making bedGraph
    awk '{ printf "%s\\t%d\\t%d\\t%2.3f\\n", \$1, \$2, \$3, \$9 }' \\
    ${prefix}.tsv \\
    > ${prefix}.bdg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.tsv
    touch  ${prefix}.bdg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
