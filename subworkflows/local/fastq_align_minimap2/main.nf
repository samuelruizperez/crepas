//
// Alignment with MINIMAP2 
//

include { MINIMAP2_ALIGN                 } from '../../../modules/nf-core/minimap2/align/main'
include { BAM_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_stats_samtools/main'

workflow FASTQ_ALIGN_MINIMAP2 {
    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), path(index) ]
    ch_fasta        // channel (optional) : [ val(meta3), path(fasta) ]
    bam_format
    bam_index_extension
    cigar_paf_format
    cigar_bam

    main:
    ch_versions = channel.empty()

    //
    // Map reads with MINIMAP2_ALIGN
    //

    MINIMAP2_ALIGN (
        ch_reads,
        ch_index,
        bam_format,
        bam_index_extension,
        cigar_paf_format,
        cigar_bam
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions_minimap2.first())

    MINIMAP2_ALIGN
        .out
        .bam
        .join(MINIMAP2_ALIGN.out.index, by: 0)
        .set { ch_bam_bai }

    //
    // Run samtools stats, flagstat and idxstats
    //
    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = MINIMAP2_ALIGN.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = MINIMAP2_ALIGN.out.index    // channel: [ val(meta), path(bai) ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
