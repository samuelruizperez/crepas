
//
// Call peaks with epic2, annotate with HOMER and perform downstream QC
//

include { BAM_REMOVE_SCAFFOLDS     } from '../../../modules/local/bam_remove_scaffolds/main'
include { DANPOS2_DPEAK           } from '../../../modules/local/danpos2/dpeak/main'


workflow BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Remove scaffolds from BAM files
    //
    BAM_REMOVE_SCAFFOLDS (
        ch_bam
    )
    ch_versions = ch_versions.mix(BAM_REMOVE_SCAFFOLDS.out.versions.first())


    // Branch channels based on if input control is present
    BAM_REMOVE_SCAFFOLDS.out.bam
        .branch { meta, bam ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, bam]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta.id, meta.antibody, meta, bam]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta, bam]
        }
        .set { ch_bam_by_type }

    // For non-downsampled files, duplicate input ipcontrols for each antibody 
    ch_bam_by_type
        .ipcontrols
        .branch { id, meta, bam ->
            dsp: meta.input_control_of_antibody && meta.dSp_total_mapped_reads
                return [id, meta.input_control_of_antibody, meta, bam]
            not_dsp: !meta.input_control_of_antibody && !meta.dSp_total_mapped_reads
                return [id, meta, bam]
        }
        .set { ch_bam_ipcontrols }
    
    ch_bam_ipcontrols
        .not_dsp
        .combine(ch_bam_by_type.ips_with_ipcontrol, by: 0) // combine by control id only
        .map { control_id, control_meta, control_bam, ip_antibody, ip_meta, ip_bam ->
            def meta_clone = control_meta.clone()
            meta_clone.input_control_of_antibody = ip_antibody
            [ control_id, meta_clone.input_control_of_antibody, meta_clone, control_bam ]
        }
        .unique()
        .set { ch_bam_ipcontrols_not_dsp }

    // Create channel: [ meta, [ip_bams_merged_reps], [ipcontrol_bams_merged_reps] ]
    ch_bam_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_ipcontrols.dsp.mix(ch_bam_ipcontrols_not_dsp), by: [0, 1])
        .map { control_id, antibody, ip_meta, ip_bam, control_meta, control_bam ->
            [ control_id, antibody, ip_meta, ip_bam, control_bam ]
        }
        .mix(ch_bam_by_type.ips_wo_ipcontrol)
        // ips_wo_ipcontrol do not have control_bam (it[4])
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            [ meta_clone.id, it[1], meta_clone, it[3], it[4] ?: [] ]
        }
        .groupTuple(by: [0, 1])
        .map {
            id, antibody, metas, ip_bams, ipcontrol_bams ->
                [ metas[0], ip_bams.flatten(), ipcontrol_bams.flatten() ]
        }
        .set { ch_ip_control_bam_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_merged_reps
        .map {
            meta, ip_bams, ipcontrol_bams ->
                "${meta}\t${ip_bams}\t${ipcontrol_bams}"
        }
        .collectFile( name: 'ch_ip_control_bam_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER" )

    ch_ip_control_bam_merged_reps
        .map {
            meta, ip_bams, ipcontrol_bams ->
                [ meta, ip_bams, ipcontrol_bams, [], [], [], [] ]

        }
        .set { ch_ip_control_bam_merged_reps }

    //
    // MODULE: Call peaks with danpos2
    //
    DANPOS2_DPEAK (
        ch_ip_control_bam_merged_reps
    )
    ch_versions = ch_versions.mix(DANPOS2_DPEAK.out.versions.first())

    emit:

    peaks                        = DANPOS2_DPEAK.out.peaks                  // channel: [ meta, "result/*.peaks.integrative.xls" ]

    versions                     = ch_versions                      // channel: [ versions.yml ]
}
