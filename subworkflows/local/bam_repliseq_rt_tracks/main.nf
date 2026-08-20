include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../../modules/nf-core/deeptools/multibamsummary/main'
include { REPLISEQ_RTNORMALIZE      } from '../../../modules/local/repliseq_rtnormalize/main'
include { FILE_SORT                 } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG     } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { FIND_CONCATENATE as REPLISEQ_RT_MULTIQC } from '../../../modules/nf-core/find/concatenate/main'

workflow BAM_REPLISEQ_RT_TRACKS {

    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai) ]
    ch_chrom_sizes_endo // channel: [ val(meta), path(chrom_sizes) ]
    ch_blacklist        // channel: [ val(meta), path(blacklist) ]
    ch_rt_header        // channel: path(repliseq_rt_header.txt)

    main:

    //
    // Group all replicates (every fraction) of the same condition together, sorted
    // deterministically by (rt_phase, brep) so bams/bais/phases/breps stay aligned.
    //
    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            // The id at this point is "<exp_type>_<sample>_<rt_phase>_bRep_<n>"
            meta_clone.id = "${meta_clone.exp_type}_${meta_clone.rt_condition}"
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
        ch_blacklist
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

    //
    // MODULE: Prepend the MultiQC custom-content header to each sample's summary row.
    // FIND_CONCATENATE concatenates its inputs in filename order rather than in the order given
    // here, so the two names have to sort the right way round: "repliseq_rt_header.txt" before
    // "rt_summary.tsv". Renaming either file without preserving that would silently put the data
    // row above the header and MultiQC would stop recognising the section.
    //
    REPLISEQ_RTNORMALIZE.out.summary
        .combine(ch_rt_header)
        .map { meta, summary, header -> [ meta, [ header, summary ] ] }
        .set { ch_rt_mqc }

    REPLISEQ_RT_MULTIQC (
        ch_rt_mqc
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

    // Only emitted when a third fraction was supplied, so these channels are empty when only early/late are inputted
    ch_rt_index_raw = REPLISEQ_RTNORMALIZE.out.index_raw.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_track_type = 'index_raw'
        [ meta_clone, bdg ]
    }

    ch_rt_index_smooth = REPLISEQ_RTNORMALIZE.out.index_smooth.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_track_type = 'index_smooth'
        [ meta_clone, bdg ]
    }

    //
    // MODULE: Lexicographically sort the RT tracks (required by UCSC_BEDGRAPHTOBIGWIG)
    //
    FILE_SORT (
        ch_rt_raw.mix(ch_rt_smooth).mix(ch_rt_index_raw).mix(ch_rt_index_smooth),
        'bedGraph'
    )

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
    mqc      = REPLISEQ_RT_MULTIQC.out.file_out // channel: [ val(meta), path(rt_mqc.tsv) ]
}
