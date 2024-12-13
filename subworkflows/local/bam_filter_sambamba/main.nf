/*
 * Filter BAM file
 */

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { SAMBAMBA_VIEW }            from '../../../modules/local/sambamba/view/main'
include { BAM_SORT_STATS_SAMTOOLS  } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'

workflow BAM_FILTER_SAMBAMBA {
    take:
    ch_bam              // channel: [ val(meta), [ bam ]]
    ch_bed                    // channel: [ bed ]
    ch_fasta                  // channel: [ fasta ]

    main:
    ch_versions = Channel.empty()

    // TODO: this one might not be necessary
    SAMTOOLS_INDEX(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMBAMBA_VIEW(
        ch_bam.join(SAMTOOLS_INDEX.out.index, by: [0]),
        ch_bed
    )

    BAM_SORT_STATS_SAMTOOLS(SAMBAMBA_VIEW.out.bam, ch_fasta)

    ch_versions = ch_versions.mix(SAMBAMBA_VIEW.out.versions,
                    BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    index    = BAM_SORT_STATS_SAMTOOLS.out.index              // channel: [ val(meta), [ index ] ]
    //bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
