include { SAMTOOLS_SORT               } from '../../../modules/nf-core/samtools/sort/main'
include { BEDTOOLS_BAMTOBED           } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { BED_TO_TAGALIGN            } from '../../../modules/local/bed_to_tagalign/main'
include { TAGALIGN_SELF_PSEUDOREPLICATES } from '../../../modules/local/tagalign_self_pseudoreplicates/main'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_fasta                          // channel: [ fasta ]
    ch_gtf                            // channel: [ gtf ]
    ch_chrom_sizes                    // channel: [ bed ]
    macs_gsize                        // integer: value for --macs_gsize parameter



    main:

    ch_versions = Channel.empty()

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
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())

    }

    //
    // MODULE: Convert BED to TAGALIGN
    //
    BED_TO_TAGALIGN (
        BEDTOOLS_BAMTOBED.out.bed
    )
    ch_versions = ch_versions.mix(BED_TO_TAGALIGN.out.versions.first())

    //
    // MODULE: Generate self-pseudoreplicates from TAGALIGN files
    //
    TAGALIGN_SELF_PSEUDOREPLICATES (
        BED_TO_TAGALIGN.out.tagalign
    )
    ch_versions = ch_versions.mix(TAGALIGN_SELF_PSEUDOREPLICATES.out.versions.first())

    // Add pseudoreplicate to metadata
    TAGALIGN_SELF_PSEUDOREPLICATES.out.tagalign
        .map { meta, tagalign ->
            def meta_clone = meta.clone()
            meta_clone.is_pseudoreplicate = true
            [ meta_clone, tagalign ]
        }
        .set {ch_self_pseudoreps}
    
    // Create
    BEDTOOLS_BAMTOBED.out.bed
        .mix(ch_self_pseudoreps)
        .branch { meta, tagalign ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, tagalign]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta.id, meta.antibody, meta, tagalign]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta, tagalign]
        }
        .set { ch_ta_by_type }

        // For non-downsampled files, duplicate input ipcontrols for each antibody 
        ch_ta_by_type
            .ipcontrols
            .branch { id, meta, tagalign ->
                dsp: meta.input_control_of_antibody && meta.dSp_total_mapped_reads
                    return [id, meta.input_control_of_antibody, meta, tagalign]
                not_dsp: !meta.input_control_of_antibody && !meta.dSp_total_mapped_reads
                    return [id, meta, tagalign]
            }
            .set { ch_ta_ipcontrols }

        ch_ta_ipcontrols
            .not_dsp
            .combine(ch_ta_by_type.ips_with_ipcontrol, by: 0) // combine by control id only
            .map { ipcontrol_id, ipcontrol_meta, ipcontrol_ta, ip_antibody, ip_meta, ip_ta ->
                def meta_clone = ipcontrol_meta.clone()
                meta_clone.input_control_of_antibody = ip_antibody
                [ ipcontrol_id, meta_clone.input_control_of_antibody, meta_clone, ipcontrol_ta ]
            }
            .unique()
            .set { ch_ta_ipcontrols_not_dsp }

    // Create channel: [ meta, [ip_tas_merged_reps], [ipcontrol_tas_merged_reps] ]
    ch_ta_by_type
        .ips_with_ipcontrol
        .combine(ch_ta_ipcontrols.dsp.mix(ch_ta_ipcontrols_not_dsp), by: [0, 1])
        .map { ipcontrol_id, antibody, ip_meta, ip_ta, ipcontrol_meta, ipcontrol_ta ->
            [ ipcontrol_id, antibody, ip_meta, ip_ta, ipcontrol_ta ]
        }
        .mix(ch_ta_by_type.ips_wo_ipcontrol)
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            [meta_clone.id, it[1], meta_clone, it[3], it[4] ?: []]
        }
        .groupTuple(by: [0, 1])
        .map { id, antibody, metas, ip_tas, ipcontrol_tas ->
            [metas[0], ip_tas.flatten(), ipcontrol_tas.flatten()]
        }
        .set { ch_tas_reps_and_pseudoreps }


    TAGALIGN_POOL (
        ch_tas_reps_and_pseudoreps
    )
    ch_versions = ch_versions.mix(TAGALIGN_POOL.out.versions.first())

    emit:
    peaks                        = ch_macs3_peaks                   // channel: [ val(meta), [ peaks ] ]


    versions                     = ch_versions                      // channel: [ versions.yml ]
}
