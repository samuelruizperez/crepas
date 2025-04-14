include { BAM_STATS_SAMTOOLS as BAM_STATS_SAMTOOLS_INITIAL } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { BAM_FLAGSTAT_MAPPED } from '../../../modules/local/bam_flagstat_mapped/main'
include { PICARD_DOWNSAMPLESAM } from '../../../modules/local/picard/downsamplesam/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS as BAM_STATS_SAMTOOLS_FINAL } from '../../../subworkflows/nf-core/bam_stats_samtools/main'

workflow BAM_DOWNSAMPLE {

    take:
    ch_bam_bai        // channel: [ val(meta), [ bam ] , [ bai ] ]
    ch_fasta          // channel: [ val(meta), [ fasta ] ]
    ch_fai            // channel: [ val(meta), [ fai ] ]
    downsampling_method // string: e.g. 'min_by_type'

    main:

    ch_versions = Channel.empty()

    BAM_STATS_SAMTOOLS_INITIAL (
        ch_bam_bai,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS_INITIAL.out.versions)

    //
    // MODULE: Add total_mapped_reads to each bam meta
    //
    BAM_FLAGSTAT_MAPPED (
        BAM_STATS_SAMTOOLS_INITIAL.out.flagstat
    )
    ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED.out.versions)

    //
    // MODULE: Extract total_mapped_reads from flagstat
    //
    BAM_FLAGSTAT_MAPPED.out.txt
        .map {
            meta, total ->
                [ meta.id, total.splitCsv(header:false)[0][0] ]
        }
        .set { ch_total }

    // Add the total_mapped_reads to the bam meta
    ch_bam_bai
        .map {
            meta, bam, bai ->
                [ meta.id, meta, bam, bai ]
        }
        .combine(ch_total, by: 0)
        .map {
            id, meta, bam, bai, total ->
                meta_clone = meta.clone()
                meta_clone.total_mapped_reads = total.toDouble()
                [ meta_clone, bam, bai ]
        }
        .set { ch_bam_bai }

    if (downsampling_method == 'min_by_type') {

        // Split ChIPs and inputs
        ch_bam_bai
            .branch { meta, bam, bai ->
                ips: !meta.is_control
                controls: meta.is_control
            }
            .set { ch_bam_bai_b }

        // Get the minimum number of reads among the ChIPs
        ch_bam_bai_b
            .ips
            .min { it[0].total_mapped_reads }
            .map { meta, bam, bai -> meta.total_mapped_reads }
            .set { ch_min_ips }

        // Get the downsample proportion for each ChIP based on the minimum
        ch_bam_bai_b
            .ips
            .combine(ch_min_ips)
            .map { meta, bam, bai, min ->
                def meta_clone = meta.clone()
                meta_clone.downsampleA_prob = min / meta_clone.total_mapped_reads
                [ meta_clone, bam, bai ]
            }
            .set { ch_bam_bai_ips_ds }

        // Get the minimum number of reads among the inputs
        ch_bam_bai_b
            .controls
            .min { it[0].total_mapped_reads }
            .set { ch_min_controls }

        // Get the downsample proportion for each input based on the minimum
        ch_bam_bai_b
            .controls
            .combine(ch_min_controls)
            .map { meta, bam, bai, min ->
                def meta_clone = meta.clone()
                meta_clone.downsampleA_prob = min / meta_clone.total_mapped_reads
                [ meta_clone, bam, bai ]
            }
            .set { ch_bam_bai_controls_ds }

        // Merge back the downsampled ChIPs and inputs
        ch_bam_bai_ips_ds
            .mix(ch_bam_bai_controls_ds)
            .set { ch_bam_bai_ds }

        // TODO: Remove this: Print channel for debugging
        ch_bam_bai_ds
            .map {
                meta, bam, bai ->
                    "${meta}\t${bam}\t${bai}"
            }
            .collectFile( name: 'ch_bam_bai_ds.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    }

    //
    // MODULE: Downsample BAMs
    //
    PICARD_DOWNSAMPLESAM (
        ch_bam_bai_ds,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(PICARD_DOWNSAMPLESAM.out.versions)

    //
    // MODULE: Index BAMs
    //
    SAMTOOLS_INDEX (
        PICARD_DOWNSAMPLESAM.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //
    // SUBWORKFLOW: Run SAMtools stats, flagstat and idxstats
    //
    BAM_STATS_SAMTOOLS_FINAL (
        PICARD_DOWNSAMPLESAM.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS_FINAL.out.versions)


    emit:
    bam      = PICARD_DOWNSAMPLESAM.out.bam   // channel: [ val(meta), [ bam ] ]
    index    = SAMTOOLS_INDEX.out.index           // channel: [ val(meta), [ index ] ]
    stats    = BAM_STATS_SAMTOOLS_FINAL.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS_FINAL.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS_FINAL.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}

