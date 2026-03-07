process SPIKEIN_BARCODE_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99b2a7149c1943265175ef013ca8c247c9847b14a0d2e802d0cbdda4a58458a5/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_r-argparse_r-ggpmisc_pruned:2c7c689513df97b4' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.tsv"), emit: merged_counts
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.merged_counts"
    if (counts.size() > 1) {
        """
        spikein_barcode_merge.R \\
            ${args} \\
            --count_tables ${counts.join(' ')} \\
            --output ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        ln -s ${counts[0]} ${prefix}.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.merged_counts"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
