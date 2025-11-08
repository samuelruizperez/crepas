include { PICARD_DOWNSAMPLESAM                           } from '../../../modules/local/picard/downsamplesam/main'
include { SAMTOOLS_INDEX                                 } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                             } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { BAM_FLAGSTAT_MAPPED as BAM_FLAGSTAT_MAPPED_DSP } from '../../../modules/local/bam_flagstat_mapped/main'

workflow BAM_DOWNSAMPLE {
    take:
    ch_bam_bai                  // channel: [ val(meta), [ bam ] , [ bai ] ]
    ch_fasta                    // channel: [ val(meta), [ fasta ] ]
    ch_fai                      // channel: [ val(meta), [ fai ] ]
    genome                      // string: e.g. 'hg38'
    spikein_genome              // string: e.g. 'mm10'
    downsampling_method         // string: e.g. 'min_by_type'
    downsampling_endo_threshold // int: e.g. 10000000
    downsampling_exo_threshold  // int: e.g. 300000

    main:

    ch_versions = channel.empty()

    // TODO: simplify this subworkflow: several steps are repeated for the different downsampling methods
    // TODO: check if samples without input control and samples without spike-in can still be downsampled
    // TODO: perhaps add:
    // 'min_total_across' downsamples all ChIPs and inputs (both endo and exo) of the same antibody to the minimum number of total reads (endo + exo) across all ChIPs and inputs of the same antibody.
    // 'min_total_by_type' downsamples all ChIPs (both endo and exo) to the minimum number of total reads (endo + exo) among the ChIPs, and all inputs to the minimum


    // 'min_endo_across' downsamples all ChIPs and inputs (both endo and exo) of the same antibody to the minimum number of endogenous reads across all ChIPs and inputs of the same antibody.
    if (downsampling_method == 'min_endo_across') {

        ch_bam_bai
            .branch { meta, bam, bai ->
                endo_ip: meta.genome == genome && !meta.is_input_control
                return [meta.input_control, meta, bam, bai]
                endo_ipcontrol: meta.genome == genome && meta.is_input_control
                return [meta.id, meta, bam, bai]
                exo_ip: meta.genome == spikein_genome && !meta.is_input_control
                return [meta.id, meta, bam, bai]
                exo_ipcontrol: meta.genome == spikein_genome && meta.is_input_control
                return [meta.id, meta, bam, bai]
            }
            .set { ch_bam_bai_genome_type }

        ch_bam_bai_genome_type.endo_ip
            .mix(ch_bam_bai_genome_type.endo_ipcontrol)
            .map { ipcontrol_id, meta, bam, bai ->
                def total = meta[meta.ref_total_mapped_reads_for_dSp]
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [meta.exp_type, antibody, total, meta, bam, bai]
            }
            .groupTuple(by: [0, 1])
            .map { exp_type, antibody, totals, metas, bams, bais ->
                // min_endo should be the minimum number in totals above the downsampling_endo_threshold
                def filtered_totals = totals.findAll { it -> it >= downsampling_endo_threshold }
                // If there are no totals above the threshold, we use the minimum of all totals
                def min_endo = filtered_totals ? filtered_totals.min() : totals.min()
                // These two lines are just to keep track of the dSp reference sample
                def min_endo_id = metas[totals.indexOf(min_endo)].id
                def min_endo_genome = metas[totals.indexOf(min_endo)].genome
                def downsampling_ref_total_key = metas[totals.indexOf(min_endo)].find { it -> it.value == min_endo }.key
                [exp_type, antibody, min_endo, min_endo_id, min_endo_genome, downsampling_ref_total_key, metas, bams, bais]
            }
            .transpose()
            .map { exp_type, antibody, min_endo, min_endo_id, min_endo_genome, downsampling_ref_total_key, meta, bam, bai ->
                def meta_clone = meta.clone()
                def downsampling_prob = min_endo / meta_clone[downsampling_ref_total_key]
                // If downsampling_prob is higher than 1 (when min_endo is higher due to downsampling_endo_threshold), we set it to 1
                meta_clone.downsampling_prob = downsampling_prob > 1 ? 1 : downsampling_prob
                meta_clone.downsampling_ref_total = min_endo
                meta_clone.downsampling_ref_sample = min_endo_id
                meta_clone.downsampling_ref_sample_genome = min_endo_genome
                meta_clone.downsampling_ref_total_key = downsampling_ref_total_key
                [meta_clone.id, meta_clone, bam, bai]
            }
            .set { ch_bam_bai_endo }

        // Now we copy the downsampling probability (by meta.id) to the exogenous bams and bais,
        // since we want the ratio of endo and exo reads in each sample to stay the same after downsampling
        ch_bam_bai_endo
            .combine(ch_bam_bai_genome_type.exo_ip.mix(ch_bam_bai_genome_type.exo_ipcontrol), by: 0)
            .map { id, endo_meta, endo_bam, endo_bai, exo_meta, exo_bam, exo_bai ->
                def meta_clone = exo_meta.clone()
                meta_clone.downsampling_prob = endo_meta.downsampling_prob
                meta_clone.downsampling_ref_total = endo_meta.downsampling_ref_total
                meta_clone.downsampling_ref_sample = endo_meta.downsampling_ref_sample
                meta_clone.downsampling_ref_sample_genome = endo_meta.downsampling_ref_sample_genome
                meta_clone.downsampling_ref_total_key = endo_meta.downsampling_ref_total_key
                // propagate input_control_of_antibody if present (i.e., in controls)
                if (endo_meta.containsKey('input_control_of_antibody')) {
                    meta_clone.input_control_of_antibody = endo_meta.input_control_of_antibody
                }
                [meta_clone.id, meta_clone, exo_bam, exo_bai]
            }
            .set { ch_bam_bai_exo }

        ch_bam_bai_endo
            .mix(ch_bam_bai_exo)
            .map { id, meta, bam, bai ->
                def meta_clone = meta.clone()
                [meta_clone + [downsampling_method: 'min_endo_across'], bam, bai]
            }
            .set { ch_bam_bai_to_ds }
    }
    else if (downsampling_method == 'min_exo_across') {

        ch_bam_bai
            .branch { meta, bam, bai ->
                endo_ip: meta.genome == genome && !meta.is_input_control
                return [meta.id, meta, bam, bai]
                endo_ipcontrol: meta.genome == genome && meta.is_input_control
                return [meta.id, meta, bam, bai]
                exo_ip: meta.genome == spikein_genome && !meta.is_input_control
                return [meta.input_control, meta, bam, bai]
                exo_ipcontrol: meta.genome == spikein_genome && meta.is_input_control
                return [meta.id, meta, bam, bai]
            }
            .set { ch_bam_bai_genome_type }

        ch_bam_bai_genome_type.exo_ip
            .mix(ch_bam_bai_genome_type.exo_ipcontrol)
            .map { ipcontrol_id, meta, bam, bai ->
                def total = meta[meta.ref_total_mapped_reads_for_dSp]
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [meta.exp_type, antibody, total, meta, bam, bai]
            }
            .groupTuple(by: [0, 1])
            .map { exp_type, antibody, totals, metas, bams, bais ->
                // min_exo should be the minimum number in totals above the downsampling_exo_threshold
                def filtered_totals = totals.findAll { it -> it >= downsampling_exo_threshold }
                def min_exo = filtered_totals ? filtered_totals.min() : totals.min()
                // These two lines are just to keep track of the dSp reference sample
                def min_exo_id = metas[totals.indexOf(min_exo)].id
                def min_exo_genome = metas[totals.indexOf(min_exo)].genome
                def downsampling_ref_total_key = metas[totals.indexOf(min_exo)].find { it -> it.value == min_exo }.key
                [exp_type, antibody, min_exo, min_exo_id, min_exo_genome, downsampling_ref_total_key, metas, bams, bais]
            }
            .transpose()
            .map { exp_type, antibody, min_exo, min_exo_id, min_exo_genome, downsampling_ref_total_key, meta, bam, bai ->
                def meta_clone = meta.clone()
                def downsampling_prob = min_exo / meta_clone[downsampling_ref_total_key]
                // If downsampling_prob is higher than 1 (when min_exo is higher due to downsampling_exo_threshold), we set it to 1
                meta_clone.downsampling_prob = downsampling_prob > 1 ? 1 : downsampling_prob
                meta_clone.downsampling_ref_total = min_exo
                meta_clone.downsampling_ref_sample = min_exo_id
                meta_clone.downsampling_ref_sample_genome = min_exo_genome
                meta_clone.downsampling_ref_total_key = downsampling_ref_total_key
                [meta_clone.id, meta_clone, bam, bai]
            }
            .set { ch_bam_bai_exo }

        // Now we copy the downsampling probability by meta.id to the exogenous bams and bais:
        // We want the ratio of endo and exo reads in each sample to stay the same after downsampling
        ch_bam_bai_exo
            .combine(ch_bam_bai_genome_type.endo_ip.mix(ch_bam_bai_genome_type.endo_ipcontrol), by: 0)
            .map { id, exo_meta, exo_bam, exo_bai, endo_meta, endo_bam, endo_bai ->
                def meta_clone = endo_meta.clone()
                meta_clone.downsampling_prob = exo_meta.downsampling_prob
                meta_clone.downsampling_ref_total = exo_meta.downsampling_ref_total
                meta_clone.downsampling_ref_sample = exo_meta.downsampling_ref_sample
                meta_clone.downsampling_ref_sample_genome = exo_meta.downsampling_ref_sample_genome
                meta_clone.downsampling_ref_total_key = exo_meta.downsampling_ref_total_key
                // propagate input_control_of_antibody if present (i.e., in controls)
                if (exo_meta.containsKey('input_control_of_antibody')) {
                    meta_clone.input_control_of_antibody = exo_meta.input_control_of_antibody
                }
                [meta_clone.id, meta_clone, endo_bam, endo_bai]
            }
        set { ch_bam_bai_endo }

        ch_bam_bai_exo
            .mix(ch_bam_bai_endo)
            .map { id, meta, bam, bai ->
                def meta_clone = meta.clone()
                [meta_clone + [downsampling_method: 'min_exo_across'], bam, bai]
            }
            .set { ch_bam_bai_to_ds }
    }
    else if (downsampling_method == 'min_endo_by_type') {

        ch_bam_bai
            .branch { meta, bam, bai ->
                endo_ip: meta.genome == genome && !meta.is_input_control
                return [meta.input_control, meta, bam, bai]
                endo_ipcontrol: meta.genome == genome && meta.is_input_control
                return [meta.id, meta, bam, bai]
                exo_ip: meta.genome == spikein_genome && !meta.is_input_control
                return [meta.id, meta, bam, bai]
                exo_ipcontrol: meta.genome == spikein_genome && meta.is_input_control
                return [meta.id, meta, bam, bai]
            }
            .set { ch_bam_bai_genome_type }

        ch_bam_bai_genome_type.endo_ip
            .mix(ch_bam_bai_genome_type.endo_ipcontrol)
            .map { ipcontrol_id, meta, bam, bai ->
                def total = meta[meta.ref_total_mapped_reads_for_dSp]
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [meta.exp_type, antibody, meta.is_input_control, total, meta, bam, bai]
            }
            .groupTuple(by: [0, 1, 2])
            .map { exp_type, antibody, is_input_control, totals, metas, bams, bais ->
                // min_endo should be the minimum number in totals above the downsampling_endo_threshold
                def filtered_totals = totals.findAll { it -> it >= downsampling_endo_threshold }
                // If there are no totals above the threshold, we use the minimum of all totals
                def min_endo = filtered_totals ? filtered_totals.min() : totals.min()
                def min_endo_id = metas[totals.indexOf(min_endo)].id
                def min_endo_genome = metas[totals.indexOf(min_endo)].genome
                // get the key name of the meta item from which we got the min_endo (e.g. "flT3_total_mapped_reads")
                def downsampling_ref_total_key = metas[totals.indexOf(min_endo)].find { it -> it.value == min_endo }.key
                [exp_type, antibody, is_input_control, min_endo, min_endo_id, min_endo_genome, downsampling_ref_total_key, metas, bams, bais]
            }
            .transpose()
            .map { exp_type, antibody, is_input_control, min_endo, min_endo_id, min_endo_genome, downsampling_ref_total_key, meta, bam, bai ->
                def meta_clone = meta.clone()
                def downsampling_prob = min_endo / meta_clone[downsampling_ref_total_key]
                // If downsampling_prob is higher than 1 (when min_endo is higher due to downsampling_endo_threshold), we set it to 1
                meta_clone.downsampling_prob = downsampling_prob > 1 ? 1 : downsampling_prob
                meta_clone.downsampling_ref_total = min_endo
                // These two lines are just to keep track of the dSp reference sample
                meta_clone.downsampling_ref_sample = min_endo_id
                meta_clone.downsampling_ref_sample_genome = min_endo_genome
                meta_clone.downsampling_ref_total_key = downsampling_ref_total_key
                [meta_clone.id, meta_clone, bam, bai]
            }
            .set { ch_bam_bai_endo }

        // Now we copy the downsampling probability by meta.id to the exogenous bams and bais:
        // We want the ratio of endo and exo reads in each sample to stay the same after downsampling
        ch_bam_bai_endo
            .combine(ch_bam_bai_genome_type.exo_ip.mix(ch_bam_bai_genome_type.exo_ipcontrol), by: 0)
            .map { id, endo_meta, endo_bam, endo_bai, exo_meta, exo_bam, exo_bai ->
                def meta_clone = exo_meta.clone()
                meta_clone.downsampling_prob = endo_meta.downsampling_prob
                meta_clone.downsampling_ref_total = endo_meta.downsampling_ref_total
                meta_clone.downsampling_ref_sample = endo_meta.downsampling_ref_sample
                meta_clone.downsampling_ref_sample_genome = endo_meta.downsampling_ref_sample_genome
                meta_clone.downsampling_ref_total_key = endo_meta.downsampling_ref_total_key
                // propagate input_control_of_antibody if present (i.e., in controls)
                if (endo_meta.containsKey('input_control_of_antibody')) {
                    meta_clone.input_control_of_antibody = endo_meta.input_control_of_antibody
                }
                [meta_clone.id, meta_clone, exo_bam, exo_bai]
            }
            .set { ch_bam_bai_exo }

        ch_bam_bai_endo
            .mix(ch_bam_bai_exo)
            .map { id, meta, bam, bai ->
                def meta_clone = meta.clone()
                [meta_clone + [downsampling_method: 'min_endo_by_type'], bam, bai]
            }
            .set { ch_bam_bai_to_ds }
    }
    else if (downsampling_method == 'min_exo_by_type') {
        ch_bam_bai
            .branch { meta, bam, bai ->
                endo_ip: meta.genome == genome && !meta.is_input_control
                return [meta.id, meta, bam, bai]
                endo_ipcontrol: meta.genome == genome && meta.is_input_control
                return [meta.id, meta, bam, bai]
                exo_ip: meta.genome == spikein_genome && !meta.is_input_control
                return [meta.input_control, meta, bam, bai]
                exo_ipcontrol: meta.genome == spikein_genome && meta.is_input_control
                return [meta.id, meta, bam, bai]
            }
            .set { ch_bam_bai_genome_type }

        ch_bam_bai_genome_type.exo_ip
            .mix(ch_bam_bai_genome_type.exo_ipcontrol)
            .map { ipcontrol_id, meta, bam, bai ->
                def total = meta[meta.ref_total_mapped_reads_for_dSp]
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [meta.exp_type, antibody, meta.is_input_control, total, meta, bam, bai]
            }
            .groupTuple(by: [0, 1, 2])
            .map { exp_type, antibody, is_input_control, totals, metas, bams, bais ->
                // min_exo should be the minimum number in totals above the downsampling_exo_threshold
                def filtered_totals = totals.findAll { it -> it >= downsampling_exo_threshold }
                // If there are no totals above the threshold, we use the minimum of all totals
                def min_exo = filtered_totals ? filtered_totals.min() : totals.min()
                // These two lines are just to keep track of the dSp reference sample
                def min_exo_id = metas[totals.indexOf(min_exo)].id
                def min_exo_genome = metas[totals.indexOf(min_exo)].genome
                def downsampling_ref_total_key = metas[totals.indexOf(min_exo)].find { it -> it.value == min_exo }.key
                [exp_type, antibody, is_input_control, min_exo, min_exo_id, min_exo_genome, downsampling_ref_total_key, metas, bams, bais]
            }
            .transpose()
            .map { exp_type, antibody, is_input_control, min_exo, min_exo_id, min_exo_genome, downsampling_ref_total_key, meta, bam, bai ->
                def meta_clone = meta.clone()
                // If downsampling_prob is higher than 1 (when min_exo is higher due to downsampling_exo_threshold), we set it to 1
                def downsampling_prob = min_exo / meta_clone[downsampling_ref_total_key]
                meta_clone.downsampling_prob = downsampling_prob > 1 ? 1 : downsampling_prob
                meta_clone.downsampling_ref_total = min_exo
                // These two lines are just to keep track of the dSp reference sample
                meta_clone.downsampling_ref_sample = min_exo_id
                meta_clone.downsampling_ref_sample_genome = min_exo_genome
                meta_clone.downsampling_ref_total_key = downsampling_ref_total_key
                [meta_clone.id, meta_clone, bam, bai]
            }
            .set { ch_bam_bai_exo }

        // Now we copy the downsampling probability by meta.id to the exogenous bams and bais:
        // We want the ratio of endo and exo reads in each sample to stay the same after downsampling
        ch_bam_bai_exo
            .combine(ch_bam_bai_genome_type.endo_ip.mix(ch_bam_bai_genome_type.endo_ipcontrol), by: 0)
            .map { id, exo_meta, exo_bam, exo_bai, endo_meta, endo_bam, endo_bai ->
                def meta_clone = endo_meta.clone()
                meta_clone.downsampling_prob = exo_meta.downsampling_prob
                meta_clone.downsampling_ref_total = exo_meta.downsampling_ref_total
                meta_clone.downsampling_ref_sample = exo_meta.downsampling_ref_sample
                meta_clone.downsampling_ref_sample_genome = exo_meta.downsampling_ref_sample_genome
                meta_clone.downsampling_ref_total_key = exo_meta.downsampling_ref_total_key
                // propagate input_control_of_antibody if present (i.e., in controls)
                if (exo_meta.containsKey('input_control_of_antibody')) {
                    meta_clone.input_control_of_antibody = exo_meta.input_control_of_antibody
                }
                [meta_clone.id, meta_clone, endo_bam, endo_bai]
            }
            .set { ch_bam_bai_endo }

        ch_bam_bai_exo
            .mix(ch_bam_bai_endo)
            .map { id, meta, bam, bai ->
                def meta_clone = meta.clone()
                [meta_clone + [downsampling_method: 'min_exo_by_type'], bam, bai]
            }
            .set { ch_bam_bai_to_ds }
    }

    // TODO: save for debugging
    ch_bam_bai_to_ds
        .map { meta, bam, bai ->
            "${meta}\t${bam}\t${bai}"
        }
        .collectFile(name: 'ch_bam_bai_to_ds.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_DOWNSAMPLE")

    //
    // MODULE: Downsample BAMs
    //
    PICARD_DOWNSAMPLESAM (
        ch_bam_bai_to_ds,
        ch_fasta,
        ch_fai
    )
    ch_ds_bam = PICARD_DOWNSAMPLESAM.out.bam
    ch_versions = ch_versions.mix(PICARD_DOWNSAMPLESAM.out.versions)

    //
    // MODULE: Index BAMs
    //
    SAMTOOLS_INDEX(
        ch_ds_bam
    )
    ch_ds_index = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //
    // SUBWORKFLOW: Run SAMtools stats, flagstat and idxstats
    //
    BAM_STATS_SAMTOOLS(
        ch_ds_bam.join(ch_ds_index, by: 0),
        ch_fasta
    )
    ch_ds_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    //
    // MODULE: Extract total mapped reads from flagstats
    //
    BAM_FLAGSTAT_MAPPED_DSP(
        ch_ds_flagstat
    )
    ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED_DSP.out.versions)

    // Extract the total mapped reads (after downsampling) from the text file
    BAM_FLAGSTAT_MAPPED_DSP.out.txt
        .map { meta, total ->
            [meta, total.splitCsv(header: false)[0][0]]
        }
        .set { ch_dsp_total }

    // Add the total_mapped_reads (after downsampling) to the bams' and bais' metas
    ch_ds_bam
        .combine(ch_ds_index, by: 0)
        .map { meta, bam, bai ->
            [meta, bam, bai]
        }
        .combine(ch_dsp_total, by: 0)
        .map { meta, bam, bai, total ->
            def meta_clone = meta.clone()
            meta_clone.dSp_total_mapped_reads = total.toDouble()
            [meta_clone, bam, bai]
        }
        .set { ch_ds_bam_bai }

    emit:
    bam      = ch_ds_bam_bai.map { meta, bam, bai -> [ meta, bam ] } // channel: [ val(meta), [ bam ] ]
    bai      = ch_ds_bam_bai.map { meta, bam, bai -> [ meta, bai ] } // channel: [ val(meta), [ index ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats                          // channel: [ val(meta), [ stats ] ]
    flagstat = ch_ds_flagstat                                        // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats                       // channel: [ val(meta), [ idxstats ] ]
    versions = ch_versions                                           // channel: [ versions.yml ]
}
