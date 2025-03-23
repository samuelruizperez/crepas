//
// Create a summary of samtools stat tables
//

include { SAMTOOLS_STATS_TRANSPOSE             } from '../../../modules/local/samtools/stats_transpose/main'
include { SAMTOOLS_STATS_CAT                  } from '../../../modules/local/samtools/stats_cat/main'

workflow SAMTOOLS_STATS_SUMMARY {

    take:
    ch_stats

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_STATS_TRANSPOSE (
        ch_stats
    )
    ch_col_stats = SAMTOOLS_STATS_TRANSPOSE.out.t_stats.collect{ it[1] }
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_TRANSPOSE.out.versions.first())

    SAMTOOLS_STATS_CAT (
        ch_col_stats
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_CAT.out.versions.first())

    emit:

    summary     = SAMTOOLS_STATS_CAT.out.cat    // channel: [ val(meta), [ summary ] ]
    
    versions    = ch_versions               // channel: [ versions.yml ]
}