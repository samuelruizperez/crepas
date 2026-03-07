//
// Extract SNAP-CUTANA™ Spike-in barcodes from CUT&RUN and CUT&Tag reactions
// More info: https://support.epicypher.com/docs/analyzing-snap-cutana-spike-in-controls-and-results-cut-and-tag
//
process FASTQ_EXTRACT_SPIKEIN_BARCODES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:22.04'
        : 'nf-core/ubuntu:22.04'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(barcode_table)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_extract"
    """
    printf "barcode_target\\tbarcode_id\\tbarcode_sequence\\tR1_count\\tR2_count\\n" > ${prefix}.tsv

    while IFS=\$'\\t' read -r barcode_target barcode_id barcode_sequence; do
        [[ "\$barcode_target" == "barcode_target" ]] && continue
        [[ -z "\$barcode_target" || -z "\$barcode_id" || -z "\$barcode_sequence" ]] && continue
        r1_count=\$(zgrep -cF "\$barcode_sequence" "${reads[0]}" || true)
        r2_count=\$(zgrep -cF "\$barcode_sequence" "${reads[1]}" || true)
        printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "\$barcode_target" "\$barcode_id" "\$barcode_sequence" "\$r1_count" "\$r2_count" >> ${prefix}.tsv
    done < "${barcode_table}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zgrep: \$(echo \$(zgrep --version | zgrep -oP '\\d+\\.\\d+' | head -1))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_extract"
    """
    touch  ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zgrep: \$(echo \$(zgrep --version | zgrep -oP '\\d+\\.\\d+' | head -1))
    END_VERSIONS
    """
}
