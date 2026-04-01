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
    count_barcodes_with_awk() {
        fastq_gz="\$1"
        barcode_tsv="\$2"
        gzip -cd "\$fastq_gz" | awk -F'\\t' -v OFS='\\t' '
            NR==FNR {
                if (FNR==1) next
                n++
                target[n] = \$1
                barcode_id[n] = \$2
                barcode_seq[n] = \$3
                next
            }
            (FNR % 4) == 2 {
                for (i=1; i<=n; i++) {
                    if (index(\$0, barcode_seq[i])) count[i]++
                }
            }
            END {
                for (i=1; i<=n; i++) {
                    print target[i], barcode_id[i], barcode_seq[i], count[i] + 0
                }
            }
        ' "\$barcode_tsv" -
    }

    count_barcodes_with_awk "${reads[0]}" "${barcode_table}" > r1_counts.tsv
    count_barcodes_with_awk "${reads[1]}" "${barcode_table}" > r2_counts.tsv

    printf "barcode_target\\tbarcode_id\\tbarcode_sequence\\tR1_count\\tR2_count\\n" > ${prefix}.tsv
    paste r1_counts.tsv r2_counts.tsv | awk -F'\\t' -v OFS='\\t' '{ print \$1, \$2, \$3, \$4, \$8 }' >> ${prefix}.tsv
    
    rm -f r1_counts.tsv r2_counts.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_extract"
    """
    touch  ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
