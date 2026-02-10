include { SAMTOOLS_SORT               } from '../../../modules/nf-core/samtools/sort/main'
include { BEDTOOLS_BAMTOBED           } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { BED_TO_TAGALIGN            } from '../../../modules/local/bed_to_tagalign/main'
include { TAGALIGN_SELF_PSEUDOREPLICATES } from '../../../modules/local/tagalign_self_pseudoreplicates/main'
include { CAT_CAT as TAGALIGN_POOL } from '../../../modules/nf-core/cat/cat/main'
include { PHANTOMPEAKQUALTOOLS } from '../../../modules/nf-core/phantompeakqualtools/main'


workflow BAM_ENCODE_PIPELINE {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_fasta                          // channel: [ fasta ]
    ctl_depth_ratio_threshold

    main:

    ch_versions = channel.empty()

    //
    // MODULE: Name-sorting BAM files
    //

    SAMTOOLS_SORT(
        ch_bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())


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

    //
    // MODULE: Generate self-pseudoreplicates of IP samples
    //
    TAGALIGN_SELF_PSEUDOREPLICATES (
        BED_TO_TAGALIGN.out.tagalign.filter { meta, tagalign -> !meta.is_input_control }
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
    
    // Create channel to pool: [ meta, tagaligns ]
    BED_TO_TAGALIGN.out.tagalign
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
            [ metas[0], tagaligns.flatten() ]
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
    ch_tagalign_pool = TAGALIGN_POOL.out.file_out.map { meta, tagalign -> [ meta + [ pooled: true ], tagalign ] }
    ch_versions = ch_versions.mix(TAGALIGN_POOL.out.versions.first())


    ch_tas_reps_and_pseudoreps
        .mix(ch_tagalign_pool)
        .branch { meta, tagalign ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, tagalign]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, tagalign]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, tagalign]
        }
        .set { ch_tagalign_by_type }

    // `ctl_depth_ratio` | 1.2 | If ratio of depth between controls is higher than this. then always use a pooled control for all replicates.
    ch_tagalign_by_type
        .ips_with_ipcontrol
        .combine(ch_tagalign_by_type.ipcontrols.filter { id, antibody, meta, tagalign -> !meta.pooled }, by: [0,1])
        .map { ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ->
            def pooled_ipcontrol_id = ipcontrol_id - ~/_bRep_.*$/
            [ pooled_ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ]
        }
        .groupTuple(by: [0,1])
        .map { pooled_ipcontrol_id, antibody, ip_metas, ip_tagaligns, ipcontrol_metas, ipcontrol_tagaligns ->
            // if depth ratio between controls is higher than ctl_depth_ratio, then use pooled control
            def ipcontrol_depths = ipcontrol_metas.collect { meta ->
                meta[meta.ref_total_mapped_reads_for_rpm]
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
            meta_clone.pooled_ipcontrol = ctl_depth_ratio_threshold_exceeded ?: false
            [ meta_clone.input_control, meta_clone.antibody, meta_clone, ip_tagalign]
        }
        .set { ch_tagalign_ips_with_ipcontrol }


    // TODO: save for debugging
    ch_tagalign_ips_with_ipcontrol
        .map { ipcontrol_id, antibody, meta, ip_tagalign ->
            "${ipcontrol_id}\t${antibody}\t${meta}\t${ip_tagalign}"
        }
        .collectFile(name: 'ch_tagalign_ips_with_ipcontrol.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")    

    // Create channel: [ meta, ip_tagalign, ipcontrol_tagalign ]
    ch_tagalign_ips_with_ipcontrol
        .combine(ch_tagalign_by_type.ipcontrols, by: [0,1])
        .map { ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ->
            [ ip_meta, ip_tagalign, ipcontrol_tagalign ]
        }
        .mix(ch_tagalign_by_type.ips_wo_ipcontrol)
        .map { it -> 
            [ it[0], it[1], it[2] ?: []]
        }
        .set { ch_ip_ipcontrol_tagalign }

    // TODO: save for debugging
    ch_ip_ipcontrol_tagalign
        .map { meta, ip_tagalign, ipcontrol_tagalign ->
            "${meta}\t${ip_tagalign}\t${ipcontrol_tagalign}"
        }
        .collectFile(name: 'ch_ip_ipcontrol_tagalign.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")    


    //
    // MODULE: Call peaks with phantompeakqualtools SPP
    //
    PHANTOMPEAKQUALTOOLS (
        ch_ip_ipcontrol_tagalign
    )
    ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions.first())


    emit:


    versions                     = ch_versions                      // channel: [ versions.yml ]
}
