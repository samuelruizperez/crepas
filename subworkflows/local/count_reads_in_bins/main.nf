
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_5KB } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_1KB } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BED_TO_SAF } from '../../../modules/local/bed_to_saf/main'
include { SUBREAD_FEATURECOUNTS  } from '../../../modules/nf-core/subread/featurecounts/main'

workflow COUNT_READS_IN_BINS {
    take:
    ch_bams             // channel: [ val(meta), [ peaks ] ]
    ch_chrom_sizes      // channel: [ val(meta), [ ip_bams ] ]

    main:

    ch_versions = Channel.empty()

    ch_windows_1kb = Channel.empty()
    BEDTOOLS_MAKEWINDOWS_1KB (
        ch_chrom_sizes
    )
    ch_windows_1kb = BEDTOOLS_MAKEWINDOWS_1KB.out.bed
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS_1KB.out.versions.first())

    ch_windows_5kb = Channel.empty()
    BEDTOOLS_MAKEWINDOWS_5KB (
        ch_chrom_sizes
    )
    ch_windows_5kb = BEDTOOLS_MAKEWINDOWS_5KB.out.bed
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS_5KB.out.versions.first())


    // Add window size to ch_windows_1kb and ch_windows_5kb metas
    ch_windows_1kb
        .map {
            meta, bed ->
                def meta_new = meta
                meta_new.window_size = '1kb'
                [ meta_new, bed ]
        }
        .set { ch_windows_1kb }

    ch_windows_5kb
        .map {
            meta, bed ->
                def meta_new = meta
                meta_new.window_size = '5kb'
                [ meta_new, bed ]
        }
        .set { ch_windows_5kb }

    ch_windows_all = ch_windows_1kb.mix(ch_windows_5kb)

    BED_TO_SAF (
        ch_windows_all
    )
    ch_windows_all = BED_TO_SAF.out.saf
    ch_versions = ch_versions.mix(BED_TO_SAF.out.versions.first())
        

    // Create channel: [ meta, bam, bed ]
    ch_bams
        .combine(ch_windows_all)
        // group by window size for featureCounts
        .groupTuple(by: [2,3])
        .map { metas_bams, bams, meta_windows, windows ->
            def meta_new = metas_bams[0] + ['window_size':meta_windows.window_size]
            meta_new.id = meta_windows.window_size
            [meta_new, bams, windows]
        }
        .set { ch_bam_windows }

    // TODO: save to file for debugging
    ch_bam_windows
        .map {
            meta, bams, windows ->
                "${meta}\t${bams}\t${windows}"
        }
        .collectFile(name: 'ch_bams_windows.txt', newLine: true, storeDir: "${params.outdir}")

    SUBREAD_FEATURECOUNTS (
        ch_bam_windows
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())


    emit:
    windows_counts  = SUBREAD_FEATURECOUNTS.out.counts      // channel: [ val(meta), [ counts ] ]
    windows_summary = SUBREAD_FEATURECOUNTS.out.summary     // channel: [ val(meta), [ summary ] ]

    versions                = ch_versions                   // channel: [ versions.yml ]
}
