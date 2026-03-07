//
// Merge barcode counts of technical replicates and create summary files for spike-in barcode QC
//

include { SPIKEIN_BARCODE_MERGE } from '../../../modules/local/spikein_barcode_merge/main'
include { SPIKEIN_BARCODE_QC    } from '../../../modules/local/spikein_barcode_qc/main'

workflow SPIKEIN_BARCODES {
    take:
    ch_barcode_count // channel: [ val(meta), path(barcode_counts.tsv) ]
    ch_uniq_total    // channel: [ val(meta), val(uniq_total) ]

    main:
    ch_versions = channel.empty()

    // Create channel: [ meta, spikein_barcode_counts, uniq_alignment_counts ]
    ch_barcode_count
        .map { meta, counts ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-3].join('_')
            def key = groupKey(meta_clone.id, meta_clone.trep_count) // trep_count defined in INPUT_CHECK subworkflow
            [key, meta_clone, counts]
        }
        .groupTuple(by: 0)
        .map { it ->
            [it[1][0], it[2].flatten()]
        }
        .set { ch_barcode_counts }


    //
    // MODULE: Merge spikein barcode counts of resequenced samples
    //
    SPIKEIN_BARCODE_MERGE (
        ch_barcode_counts
    )
    ch_versions = ch_versions.mix(SPIKEIN_BARCODE_MERGE.out.versions.first())


    ch_uniq_total = ch_uniq_total.map { meta, uniq_total -> [ meta.id, meta, uniq_total ]}

    SPIKEIN_BARCODE_MERGE
        .out
        .merged_counts
        .map { meta, counts -> [ meta.id, meta, counts ]}
        .combine(ch_uniq_total, by: 0)
        .map { id, meta, counts, uniq_total ->
            [ "group_all", meta, counts, uniq_total ]
        }
        .groupTuple()
        .map { group_all, metas, counts, uniq_totals ->
            [ metas[0], counts, uniq_totals ]
        }
        .set { ch_barcode_counts_uniq }

    //
    // MODULE: Calculate spike-in barcode QC metrics and plot
    //
    SPIKEIN_BARCODE_QC (
        ch_barcode_counts_uniq
    )
    ch_versions = ch_versions.mix(SPIKEIN_BARCODE_QC.out.versions.first())

  
    emit:
    spikein_barcode_qc = SPIKEIN_BARCODE_QC.out.tsv // channel: [ val(meta), path(tsv) ]

    versions = ch_versions // channel: [ versions.yml ]
}
