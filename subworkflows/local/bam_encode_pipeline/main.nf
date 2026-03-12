include { SAMTOOLS_SORT                                     } from '../../../modules/nf-core/samtools/sort/main'
include { BEDTOOLS_BAMTOBED                                 } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { BED_TO_TAGALIGN                                   } from '../../../modules/local/bed_to_tagalign/main'
include { TAGALIGN_SELF_PSEUDOREPLICATES                    } from '../../../modules/local/tagalign_self_pseudoreplicates/main'
include { CAT_CAT as TAGALIGN_POOL                          } from '../../../modules/nf-core/cat/cat/main'
include { PHANTOMPEAKQUALTOOLS as PHANTOMPEAKQUALTOOLS_SPP  } from '../../../modules/nf-core/phantompeakqualtools/main'
include { IDR                                               } from '../../../modules/nf-core/idr/main'

workflow BAM_ENCODE_PIPELINE {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_fasta                          // channel: [ fasta ]
    ctl_depth_ratio_threshold
    peak_type

    main:

    ch_versions = channel.empty()

    //
    // MODULE: Name-sorting BAM files
    //
    SAMTOOLS_SORT (
        ch_bam,
        ch_fasta,
        ''
    )

    // 2a section of the ENCODE 3 ChIP-seq pipeline:

    //
    // MODULE: Convert BAM to BED
    //
    BEDTOOLS_BAMTOBED (
        SAMTOOLS_SORT.out.bam
    )

    //
    // MODULE: Convert BED to TAGALIGN
    //
    BED_TO_TAGALIGN (
        BEDTOOLS_BAMTOBED.out.bed
    )
    ch_versions = ch_versions.mix(BED_TO_TAGALIGN.out.versions.first())


    BED_TO_TAGALIGN
        .out
        .tagalign
        .filter { meta, tagalign -> !meta.is_input_control } // We do not generate pseudoreplicates for input controls
        .map { meta, tagalign -> [ meta + [ is_pseudoreplicate: true ], tagalign ]}
        .set { ch_tagalign_ips_for_pseudoreps }

    //
    // MODULE: Generate self-pseudoreplicates of IP samples
    //
    TAGALIGN_SELF_PSEUDOREPLICATES (
        ch_tagalign_ips_for_pseudoreps
    )
    ch_versions = ch_versions.mix(TAGALIGN_SELF_PSEUDOREPLICATES.out.versions.first())

    // Add pseudoreplicate to metadata
    TAGALIGN_SELF_PSEUDOREPLICATES.out.tagalign1
        .map { meta, tagalign -> [ meta + [ pseudoreplicate: '1' ], tagalign ]}
        .mix(
            TAGALIGN_SELF_PSEUDOREPLICATES.out.tagalign2
                .map { meta, tagalign -> [ meta + [ pseudoreplicate: '2' ], tagalign ] }
        )
        .set {ch_self_pseudoreps}
    

    // Create channel: [ meta, tagaligns ] to pool replicates and pseudoreplicates
    BED_TO_TAGALIGN
        .out
        .tagalign
        .mix(ch_self_pseudoreps)
        .set { ch_tas_reps_and_pseudoreps }
        
    ch_tas_reps_and_pseudoreps
        .map { meta, tagalign ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            def antibody = meta.input_control_of_antibody ?: meta.antibody
            def pseudoreplicate = meta.pseudoreplicate ?: false
            [  meta_clone.id, antibody, pseudoreplicate, meta_clone, tagalign ]
        }
        .groupTuple(by: [0, 1, 2])
        .map { id, antibody, pseudoreplicate, metas, tagaligns ->
            // Sort metas and tagaligns to ensure consistent ordering for caching and resuming
            def sorted_tagaligns = tagaligns.sort { it -> it.name }
            def sorted_metas = metas.sort { meta -> meta.brep }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.is_pooled = true
            [ meta_clone, sorted_tagaligns ]
        }
        .set { ch_tas_reps_and_pseudoreps_to_pool }

    // TODO: save for debugging
    ch_tas_reps_and_pseudoreps_to_pool
        .map { meta, tagaligns ->
            "${meta}\t${tagaligns}"
        }
        .collectFile(name: 'ch_tas_reps_and_pseudoreps_to_pool.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")    

    //
    // MODULE: Pool replicates and pseudoreplicates with cat
    //
    TAGALIGN_POOL (
        ch_tas_reps_and_pseudoreps_to_pool
    )
    ch_tas_reps_and_pseudoreps_pooled = TAGALIGN_POOL.out.file_out

    // TODO: save for debugging
    ch_tas_reps_and_pseudoreps_pooled
        .map { meta, tagaligns ->
            "${meta}\t${tagaligns}"
        }
        .collectFile(name: 'ch_tas_reps_and_pseudoreps_pooled.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")    

    //
    // Create channel: [ meta, tagalign ] with metadata indicating whether to use pooled control or not for each sample
    //
    ch_tas_reps_and_pseudoreps
        .branch { meta, tagalign ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, tagalign]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, tagalign]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, tagalign]
        }
        .set { ch_tas_reps_and_pseudoreps_by_type }

    ch_tas_reps_and_pseudoreps_by_type
        .ips_with_ipcontrol
        .combine(ch_tas_reps_and_pseudoreps_by_type.ipcontrols, by: [0, 1]) // combine by ipcontrol_id and antibody
        .map { ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ->
            def pooled_ipcontrol_id = ipcontrol_id - ~/_bRep_.*$/
            [ pooled_ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ]
        }
        .groupTuple(by: [0, 1]) // group samples
        .map { pooled_ipcontrol_id, antibody, ip_metas, ip_tagaligns, ipcontrol_metas, ipcontrol_tagaligns ->
            // if depth ratio between controls is higher than ctl_depth_ratio, then use pooled control
            def ipcontrol_depths = ipcontrol_metas.collect { meta ->
                meta[meta.ref_total_mapped_reads_for_rpm_key]
            }
            def ctl_depth_max = ipcontrol_depths.max()
            def ctl_depth_min = ipcontrol_depths.min()
            def ctl_depth_ratio = ctl_depth_max / ctl_depth_min
            def ctl_depth_ratio_threshold_exceeded = ctl_depth_ratio > ctl_depth_ratio_threshold
            [ pooled_ipcontrol_id, antibody, ctl_depth_max, ctl_depth_min, ctl_depth_ratio, ctl_depth_ratio_threshold_exceeded, ip_metas, ip_tagaligns, ipcontrol_metas, ipcontrol_tagaligns ]
        }
        .transpose()
        .map { pooled_ipcontrol_id, antibody, ctl_depth_max, ctl_depth_min, ctl_depth_ratio, ctl_depth_ratio_threshold_exceeded, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ->
            def meta_clone = ip_meta.clone()
            if (ctl_depth_ratio_threshold_exceeded) {
                meta_clone.input_control = pooled_ipcontrol_id
            } 
            meta_clone.ctl_depth_max = ctl_depth_max
            meta_clone.ctl_depth_min = ctl_depth_min
            meta_clone.ctl_depth_ratio = ctl_depth_ratio
            meta_clone.ctl_depth_ratio_threshold_exceeded = ctl_depth_ratio_threshold_exceeded
            meta_clone.ipcontrol_is_pooled = ctl_depth_ratio_threshold_exceeded ?: false
            [ meta_clone, ip_tagalign]
        }
        .set { ch_tas_reps_and_pseudoreps_ips_with_ipcontrol }
    
    // TODO: save for debugging
    ch_tas_reps_and_pseudoreps_ips_with_ipcontrol
        .map { meta, ip_tagalign ->
            "${meta}\t${ip_tagalign}"
        }
        .collectFile(name: 'ch_tas_reps_and_pseudoreps_ips_with_ipcontrol.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")    

    // We remove id and antibody (used when branching above) to mix below
    ch_tas_reps_and_pseudoreps_by_type
        .ipcontrols
        .map { id, antibody, meta, tagalign ->
            [ meta, tagalign ]
        }
        .set { ch_tas_reps_and_pseudoreps_ipcontrols }
        
    // We mix back the rest with the ips now with updated pooled/non-pooled control metadata 
    ch_tas_reps_and_pseudoreps_by_type.ips_wo_ipcontrol
        .mix(ch_tas_reps_and_pseudoreps_pooled)
        .mix(ch_tas_reps_and_pseudoreps_ips_with_ipcontrol)
        .mix(ch_tas_reps_and_pseudoreps_ipcontrols)
        .set { ch_tagalign }

    // Create channel: [ meta, ip_tagalign, ipcontrol_tagalign ]
    ch_tagalign
        .branch { meta, tagalign ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, tagalign]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, tagalign]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, tagalign]
        }
        .set { ch_tagalign_by_type }

    ch_tagalign_by_type
        .ips_with_ipcontrol
        .combine(ch_tagalign_by_type.ipcontrols, by: [0, 1])
        .map { ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ->
            [ ip_meta, ip_tagalign, ipcontrol_tagalign ]
        }
        .mix(ch_tagalign_by_type.ips_wo_ipcontrol)
        .map { it -> [ it[0], it[1], it[2] ?: []] }
        .set { ch_tagalign_for_spp }

    // TODO: save for debugging
    ch_tagalign_for_spp
        .map { meta, ip_tagalign, ipcontrol_tagalign ->
            "${meta}\t${ip_tagalign}\t${ipcontrol_tagalign}"
        }
        .collectFile(name: 'ch_tagalign_for_spp.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")    

    //
    // MODULE: Call peaks with phantompeakqualtools SPP
    //
    PHANTOMPEAKQUALTOOLS_SPP (
        ch_tagalign_for_spp
    )
    ch_spp_peaks = PHANTOMPEAKQUALTOOLS_SPP.out.regionpeak
    ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS_SPP.out.versions.first())


    // Create channel: [ meta, [peaks1, peaks2], pooled_peaks ]
    ch_spp_peaks
        .branch { meta, peak ->
            true_replicates: !meta.pseudoreplicate && !meta.is_pooled
                return [ meta.id, meta, peak ]
            pooled_replicates: meta.is_pooled && !meta.pseudoreplicate
                return [ meta.id, peak ]
            pseudoreplicates: meta.pseudoreplicate && !meta.is_pooled
                return [ meta.id, meta, peak ]
            pooled_pseudoreplicates: meta.pseudoreplicate && meta.is_pooled
                return [ meta.id, meta, peak ]
        }
        .set { ch_spp_peaks_by_type }

    // True replicates
    ch_spp_peaks_by_type
        .true_replicates
        .map { id, meta, peak ->
            def pooled_id = id - ~/_bRep_.*$/
            [ pooled_id, meta, peak ]
        }
        .set { ch_spp_peaks_true_reps }

    ch_spp_peaks_true_reps
        .combine(ch_spp_peaks_true_reps, by: 0)
        // Ensure unordered combinations without repetition
        // See https://github.com/nextflow-io/nextflow/discussions/2109#discussioncomment-12501996
        .filter { id, meta1, peak1, meta2, peak2 -> meta1.id < meta2.id }
        .map{ id, meta1, peak1, meta2, peak2 ->
            def meta_clone = meta1.clone()
            meta_clone.id = id
            meta_clone.idr_pair_type = 'true_replicate'
            [ meta_clone.id, meta_clone, [ peak1, peak2 ] ]
        }
        .combine(ch_spp_peaks_by_type.pooled_replicates, by: 0)
        .map { id, meta, peaks, pooled_peak ->
            [ meta, peaks, pooled_peak ]
        }
        .set { ch_spp_peaks_true_reps_for_idr }

    // TODO: save for debugging
    ch_spp_peaks_true_reps_for_idr
        .map { meta, peaks, pooled_peak ->
            "${meta}\t${peaks}\t${pooled_peak}"
        }
        .collectFile(name: 'ch_spp_peaks_true_reps_for_idr.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    // Pooled pseudoreplicates
    ch_spp_peaks_by_type
        .pooled_pseudoreplicates
        .groupTuple(by: 0)
        .combine(ch_spp_peaks_by_type.pooled_replicates, by: 0)
        .map { id, metas, peaks, pooled_peak ->
            // Sort metas and peaks to ensure consistent ordering for caching and resuming
            def sorted_peaks = peaks.sort { peak -> peak.name }
            def sorted_metas = metas.sort { meta -> meta.brep }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.id = id
            meta_clone.idr_pair_type = 'pooled_pseudoreplicate'
            [ meta_clone, sorted_peaks, pooled_peak ]
        }
        .set { ch_spp_peaks_pooled_pseudoreps_for_idr }

    // TODO: save for debugging
    ch_spp_peaks_pooled_pseudoreps_for_idr
        .map { meta, peaks, pooled_peak ->
            "${meta}\t${peaks}\t${pooled_peak}"
        }
        .collectFile(name: 'ch_spp_peaks_pooled_pseudoreps_for_idr.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    // Self-pseudoreplicates
    ch_spp_peaks_by_type
        .pseudoreplicates
        .groupTuple(by: 0)
        // id, metas, peaks
        .combine(ch_spp_peaks_by_type.true_replicates, by: 0)
        .map { id, metas, peaks, true_rep_meta, true_rep_peak ->
            // Sort metas and peaks to ensure consistent ordering for caching and resuming
            def sorted_peaks = peaks.sort { peak -> peak.name }
            def sorted_metas = metas.sort { meta -> meta.pseudoreplicate }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.id = id
            meta_clone.idr_pair_type = 'self_pseudoreplicate'
            [ meta_clone, sorted_peaks, true_rep_peak ]
        }
        .set { ch_spp_peaks_self_pseudoreps_for_idr }

    // TODO: save for debugging
    ch_spp_peaks_self_pseudoreps_for_idr
        .map { meta, peaks, true_rep_peak ->
            "${meta}\t${peaks}\t${true_rep_peak}"
        }
        .collectFile(name: 'ch_spp_peaks_self_pseudoreps_for_idr.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

        
    // Use IDR to compare all pairs of matched replicates
    // (1) True replicates narrowPeak files: ${REP1_PEAK_FILE} vs. ${REP2_PEAK_FILE} IDR results transferred to Pooled-replicates narrowPeak file  ${POOLED_PEAK_FILE}
    // (2) Pooled-pseudoreplicates: ${PPR1_PEAK_FILE} vs. ${PPR2_PEAK_FILE} IDR results transferred to Pooled-replicates narrowPeak file ${POOLED_PEAK_FILE}
    // (3) Rep1 self-pseudoreplicates: ${REP1_PR1_PEAK_FILE} vs. ${REP1_PR2_PEAK_FILE} IDR results transferred to Rep1 narrowPeak file ${REP1_PEAK_FILE}
    // (4) Rep2 self-pseudoreplicates: ${REP2_PR1_PEAK_FILE} vs. ${REP2_PR2_PEAK_FILE} IDR results transferred to Rep2 narrowPeak file ${REP2_PEAK_FILE}
    ch_spp_peaks_true_reps_for_idr
        .mix(ch_spp_peaks_pooled_pseudoreps_for_idr)
        .mix(ch_spp_peaks_self_pseudoreps_for_idr)
        .set { ch_for_idr }

        
    // TODO: save for debugging
    ch_for_idr
        .map { meta, peaks, pooled_peak ->
            "${meta}\t${peaks}\t${pooled_peak}"
        }
        .collectFile(name: 'ch_for_idr.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    //
    // MODULE: IDR analysis
    //
    IDR (
        ch_for_idr,
        peak_type
    )
    ch_versions = ch_versions.mix(IDR.out.versions.first())


    emit:

    versions                     = ch_versions                      // channel: [ versions.yml ]
}
