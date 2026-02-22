process PLOT_HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-ggplot2_pruned:9dfe088cdf18888b':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-ggplot2_pruned:d2dbf7b09452d5e3' }"

    input:
    tuple val(meta), path(annos)
    path mqc_header
    val suffix

    output:
    tuple val(meta), path("*.txt")       , emit: txt
    tuple val(meta), path("*.pdf")       , emit: pdf
    tuple val(meta), path("*.tsv")       , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_homer_annotatepeaks.r \\
        -i ${annos.join(',')} \\
        -s ${annos.join(',').replaceAll("${suffix}","")} \\
        -p $prefix \\
        $args

    find ./ -type f -name "*summary.txt" -exec cat {} \\; | cat $mqc_header - > ${prefix}.summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
