//
// Filter BAM file with SAMBAMBA, optionally remove orphan reads, and generate stats with SAMTOOLS
//

include { SAMBAMBA_VIEW                                             } from '../../../modules/local/sambamba/view/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_FLT    } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_RMO    } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'
include { SAMTOOLS_SORT as SAMTOOLS_NSORT                           } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_REMOVE_ORPHANS                                        } from '../../../modules/local/bam_remove_orphans/main'
include { BAM_FLAGSTAT_MAPPED                                       } from '../../../modules/local/bam_flagstat_mapped/main'

workflow BAM_FILTER_SAMBAMBA_RMO_STATS {
    take:
    ch_bam_bai              // channel: [ val(meta), [ bam ], [ bai ]]
    ch_bed                  // channel: [ val(meta2), [ bed ] ]
    ch_fasta                // channel: [ fasta ]
    skip_orphan_removal     // boolean
    total_mapped_reads_key  // string

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: Filter BAM file with SAMBAMBA
    //
    SAMBAMBA_VIEW (
        ch_bam_bai,
        ch_bed
    )
    ch_versions = ch_versions.mix(SAMBAMBA_VIEW.out.versions)

    //
    // MODULE: Sort BAM file and generate stats with SAMTOOLS
    //
    BAM_SORT_STATS_SAMTOOLS_FLT (
        SAMBAMBA_VIEW.out.bam,
        ch_fasta
    )
    ch_bam           = BAM_SORT_STATS_SAMTOOLS_FLT.out.bam
    ch_bai           = BAM_SORT_STATS_SAMTOOLS_FLT.out.bai
    ch_flagstat      = BAM_SORT_STATS_SAMTOOLS_FLT.out.flagstat
    ch_stats         = BAM_SORT_STATS_SAMTOOLS_FLT.out.stats
    ch_idxstats      = BAM_SORT_STATS_SAMTOOLS_FLT.out.idxstats
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS_FLT.out.stats.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS_FLT.out.flagstat.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS_FLT.out.idxstats.collect { it -> it[1] })

    if (!skip_orphan_removal) {

        // Separate single-end and paired-end files (SE do not have orphans)
        ch_bam
            .branch { meta, bam ->
                se: meta.single_end
                pe: !meta.single_end
            }
            .set { ch_bam }

        ch_bai
            .branch { meta, bai ->
                se: meta.single_end
                pe: !meta.single_end
            }
            .set { ch_bai }

        ch_flagstat
            .branch { meta, flagstat ->
                se: meta.single_end
                pe: !meta.single_end
            }
            .set { ch_flagstat }

        ch_stats
            .branch { meta, stats ->
                se: meta.single_end
                pe: !meta.single_end
            }
            .set { ch_stats }

        ch_idxstats
            .branch { meta, idxstats ->
                se: meta.single_end
                pe: !meta.single_end
            }
            .set { ch_idxstats }

        //
        // MODULE: Sort BAM files by query name
        //
        SAMTOOLS_NSORT (
            ch_bam.pe,
            ch_fasta,
            ''
        )

        //
        // MODULE: Remove orphan reads left by SAMBAMBA's filtering
        //
        BAM_REMOVE_ORPHANS (
            SAMTOOLS_NSORT.out.bam,
            true
        )
        ch_versions = ch_versions.mix(BAM_REMOVE_ORPHANS.out.versions.first())

        //
        // MODULE: Sort BAM by coordinate and generate stats
        //
        BAM_SORT_STATS_SAMTOOLS_RMO (
            BAM_REMOVE_ORPHANS.out.bam,
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS_RMO.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS_RMO.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS_RMO.out.idxstats.collect { it -> it[1] })

        // Mixing SE and PE (removed orphans) files
        ch_bam      = ch_bam.se.mix(BAM_SORT_STATS_SAMTOOLS_RMO.out.bam)
        ch_bai    = ch_bai.se.mix(BAM_SORT_STATS_SAMTOOLS_RMO.out.bai)
        ch_flagstat = ch_flagstat.se.mix(BAM_SORT_STATS_SAMTOOLS_RMO.out.flagstat)
        ch_stats    = ch_stats.se.mix(BAM_SORT_STATS_SAMTOOLS_RMO.out.stats)
        ch_idxstats = ch_idxstats.se.mix(BAM_SORT_STATS_SAMTOOLS_RMO.out.idxstats)

    }

    //
    // MODULE: Extract total mapped reads from flagstats
    //
    BAM_FLAGSTAT_MAPPED (
        ch_flagstat
    )
    ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED.out.versions)

    // Extract the total mapped reads from the text file
    BAM_FLAGSTAT_MAPPED.out.txt
        .map { meta, total ->
            [meta, total.splitCsv(header: false)[0][0]]
        }
        .set { ch_total_reads }


    // Add the total_mapped_reads to the bams' and bais' metas
    ch_bam
        .join(ch_bai, by: 0)
        .combine(ch_total_reads, by: 0)
        .map { meta, bam, bai, total ->
            def meta_clone = meta.clone()
            meta_clone[total_mapped_reads_key] = total.toDouble()
            meta_clone.ref_total_mapped_reads_key = total_mapped_reads_key
            [meta_clone, bam, bai]
        }
        .multiMap { meta, bam, bai ->
            bam: [ meta, bam ]
            bai: [ meta, bai ]
        }
        .set { ch_bam_bai }

    emit:
    bam             = ch_bam_bai.bam    // channel: [ val(meta), [ bam ] ]
    bai             = ch_bam_bai.bai    // channel: [ val(meta), [ bai ] ]
    stats           = ch_stats          // channel: [ val(meta), [ stats ] ]
    flagstat        = ch_flagstat       // channel: [ val(meta), [ flagstat ] ]
    idxstats        = ch_idxstats       // channel: [ val(meta), [ idxstats ] ]
    multiqc_files   = ch_multiqc_files  // channel: [ val(meta), [ multiqc_files ] ]
    total_reads     = ch_total_reads    // channel: [ val(meta), total_mapped_reads ]

    versions = ch_versions              // channel: [ versions.yml ]
}
