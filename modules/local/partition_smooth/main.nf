process PARTITION_SMOOTH {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::perl=5.26.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2':
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(f_tab), path(r_tab)
    val radius
    val dradius
    val zradius

    output:
    tuple val(meta), path("*.txt"), emit: rfd
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in bin
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    partition_smooth.pl \\
        $f_tab \\
        $r_tab \\
        $radius \\
        $dradius \\
        $zradius \\
        > ${prefix}.RFD.txt &

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.RFD.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}
