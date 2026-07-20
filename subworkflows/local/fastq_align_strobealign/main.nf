//
// Alignment with strobealign
//

include { STROBEALIGN           } from "../../../modules/nf-core/strobealign/main"
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../../nf-core/bam_stats_samtools/main'

workflow FASTQ_ALIGN_STROBEALIGN {
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_strobealign_index
    ch_fasta_fai      // channel: [ val(meta), path(fasta), path(fai) ]
    sort_bam          // val

    main:

    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [ meta, fasta ] }

    //
    // MODULE: Map reads with strobealign
    //
    STROBEALIGN  (
        ch_reads,
        ch_fasta,
        ch_strobealign_index,
        sort_bam
    )
    ch_bam = STROBEALIGN.out.bam

    //
    // MODULE: Index BAM file with samtools
    //
    SAMTOOLS_INDEX ( ch_bam )

    ch_bam_index = ch_bam.join(SAMTOOLS_INDEX.out.index, by: 0)


    //
    // MODULE: Run samtools stats, flagstat and idxstats
    //
    BAM_STATS_SAMTOOLS ( ch_bam_index, ch_fasta_fai )

    emit:
    bam              = ch_bam      // channel: [ val(meta), [ bam ] ]
    index            = SAMTOOLS_INDEX.out.index    // channel: [ val(meta), [ index ] ]
    stats            = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
