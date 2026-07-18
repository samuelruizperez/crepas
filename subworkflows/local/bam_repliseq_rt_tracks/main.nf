include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../../modules/nf-core/deeptools/multibamsummary/main'
include { REPLISEQ_RTNORMALIZE      } from '../../../modules/local/repliseq_rtnormalize/main'
include { FILE_SORT                 } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG     } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BAM_REPLISEQ_RT_TRACKS {

    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai) ] filtered/deduped BAMs for exp_type == 'Repli-seq', meta.rt_phase == 'early'|'late'
    ch_chrom_sizes_endo // channel: [ val(meta), path(chrom_sizes) ]
    ch_blacklist        // channel: [ val(meta), path(blacklist) ]

    main:

    ch_versions = channel.empty()

    //
    // Group all replicates (both phases) of the same condition together, sorted
    // deterministically by (rt_phase, brep) so bams/bais/phases/breps stay aligned.
    //
    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            [ meta_clone.id, meta_clone, bam, bai ]
        }
        .groupTuple(by: 0)
        .map { id, metas, bams, bais ->
            def replicates = [metas, bams, bais].transpose().sort { it -> "${it[0].rt_phase}_${it[0].brep}" }
            def sorted_metas = replicates.collect { it -> it[0] }
            def meta_clone = sorted_metas[0].clone()
            def sorted_bams = replicates.collect { it -> it[1] }
            def sorted_bais = replicates.collect { it -> it[2] }
            def sorted_phases = sorted_metas.collect { meta -> meta.rt_phase }
            def sorted_breps = sorted_metas.collect { meta -> meta.brep }
            meta_clone.remove('brep')
            meta_clone.remove('rt_phase')
            [ meta_clone, sorted_bams, sorted_bais, sorted_phases, sorted_breps ]
        }
        .set { ch_condition_bam_bai }

    ch_condition_bam_bai
        .map { meta, bams, bais, phases, breps ->
            def labels = [phases, breps].transpose().collect { phase, brep -> "${phase}_${brep}" }
            [ meta, bams, bais, labels ]
        }
        .set { ch_for_multibamsummary }

    //
    // MODULE: Count reads per genomic window (bins mode)
    // across all replicates of a condition
    //
    DEEPTOOLS_MULTIBAMSUMMARY (
        ch_for_multibamsummary,
        ch_blacklist.ifEmpty([[:], []])
    )

    ch_phases_breps = ch_condition_bam_bai.map { meta, bams, bais, phases, breps -> [ meta, phases, breps ] }


    DEEPTOOLS_MULTIBAMSUMMARY.out.raw_counts
        .join(ch_phases_breps, by: 0)
        .map { meta, counts, phases, breps -> 
            [ meta, counts, phases, breps ] 
        }
        .set { ch_phases_breps_for_rtnormalize }

    //
    // MODULE: Calculate and normalize the RT track: CPM, paired/pooled/auto replicate
    // combination, optional quantile normalization between replicates, loess/roll smoothing
    //
    REPLISEQ_RTNORMALIZE (
        ch_phases_breps_for_rtnormalize
    )

    ch_rt_raw = REPLISEQ_RTNORMALIZE.out.raw.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_track_type = 'raw'
        [ meta_clone, bdg ]
    }

    ch_rt_smooth = REPLISEQ_RTNORMALIZE.out.smooth.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_track_type = 'smooth'
        [ meta_clone, bdg ]
    }

    //
    // MODULE: Lexicographically sort the RT tracks (required by UCSC_BEDGRAPHTOBIGWIG)
    //
    FILE_SORT (
        ch_rt_raw.mix(ch_rt_smooth),
        'bedGraph'
    )
    ch_versions = ch_versions.mix(FILE_SORT.out.versions.first())

    //
    // MODULE: Convert RT tracks (raw + smoothed) to bigWig
    //
    UCSC_BEDGRAPHTOBIGWIG (
        FILE_SORT.out.sorted,
        ch_chrom_sizes_endo.map { it -> it[1] }
    )

    emit:
    bedgraph = FILE_SORT.out.sorted             // channel: [ val(meta), path(bedgraph) ]
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), path(bigwig) ]
    qc       = REPLISEQ_RTNORMALIZE.out.qc      // channel: [ val(meta), path(qc.txt) ]
    versions = ch_versions                      // channel: [ path(versions.yml) ]
}
