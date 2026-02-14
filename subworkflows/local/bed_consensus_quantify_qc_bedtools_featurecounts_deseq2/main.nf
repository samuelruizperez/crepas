
//
// Call consensus peaks with BEDTools and custom scripts, annotate with HOMER, quantify with featureCounts and QC with DESeq2
//

include { HOMER_ANNOTATEPEAKS    } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { SUBREAD_FEATURECOUNTS  } from '../../../modules/nf-core/subread/featurecounts/main'

include { MACS3_CONSENSUS        } from '../../../modules/local/macs3_consensus/main'
include { ANNOTATE_BOOLEAN_PEAKS } from '../../../modules/local/annotate_boolean_peaks/main'
include { DESEQ2_QC              } from '../../../modules/local/deseq2_qc/main'

include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS } from '../../../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE as DEEPTOOLS_PLOTPROFILE_PEAKS } from '../../../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS } from '../../../modules/nf-core/deeptools/plotheatmap/main'


workflow BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 {
    take:
    ch_peaks                            // channel: [ val(meta), [ peaks ] ]
    ch_bams                             // channel: [ val(meta), [ ip_bams ] ]
    ch_bigwigs                          // channel: [ val(meta), [ bigwigs ] ]
    ch_fasta                            // channel: [ fasta ]
    ch_gtf                              // channel: [ gtf ]
    ch_deseq2_pca_header_multiqc        // channel: [ header_file ]
    ch_deseq2_clustering_header_multiqc // channel: [ header_file ]
    is_narrow_peak                      // boolean: true/false
    skip_peak_annotation                // boolean: true/false
    skip_deseq2_qc                      // boolean: true/false
    skip_consensus_plotprofile          // boolean: true/false

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //TODO: print to fiule for debugging
    ch_peaks
        .map {
            meta, peaks ->
                "${meta}\t${peaks}"
        }
        .collectFile( name: 'ch_peaks.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2" )


    // Create channels: [ meta , [ peaks ] ]
    // Where meta = [ id:antibody, exp_type, multiple_groups:true/false, replicates_exist:true/false ]
    ch_peaks
        .map {
            meta, peak ->
                [ meta.antibody, meta.exp_type, meta.id - ~/_bRep_.*$/, meta, peak ]
        }
        .tap { ch_antibody_peaks0 }
        .groupTuple(by: [0, 1])
            .map {
                antibody, exp_type, groups, metas, peaks ->
                [
                    antibody,
                    exp_type,
                    metas,
                    groups.groupBy().collectEntries { [(it.key) : it.value.size()] },
                    peaks
                ]
            }
            .tap { ch_antibody_peaks1 }
            .map {
                antibody, exp_type, metas, groups, peaks ->
                def meta_new = metas[0].clone()
                // Set meta_new.id based on exp_type and antibody presence
                if (antibody) {
                    if (exp_type == 'ATAC-seq') {
                    meta_new.id = exp_type
                    } else {
                    meta_new.id = exp_type + '_' + antibody
                    }
                } else {
                    meta_new.id = exp_type + '_no_antibody_'
                }
                meta_new.multiple_groups = groups.size() > 1
                meta_new.replicates_exist = groups.max { it.value }.value > 1
                meta_new.aligner = aligner
                [ meta_new, peaks ]
        }
        .set { ch_antibody_peaks }


    //TODO: print to fiule for debugging
    ch_antibody_peaks
        .map {
            it ->
                "${it}"
        }
        .collectFile( name: 'antibody_peaks.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2" )

    //TODO: print to fiule for debugging
    ch_antibody_peaks1
        .map {
            it ->
                "${it}"
        }
        .collectFile( name: 'ch_antibody_peaks0.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2" )

    //TODO: print to fiule for debugging
    ch_antibody_peaks0
        .map {
            it ->
                "${it}"
        }
        .collectFile( name: 'ch_antibody_peaks1.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2" )



    //
    // Generate consensus peaks across samples
    //
    MACS3_CONSENSUS (
        ch_antibody_peaks,
        is_narrow_peak
    )
    ch_versions = ch_versions.mix(MACS3_CONSENSUS.out.versions)

    //
    // Annotate consensus peaks
    //
    if (!skip_peak_annotation) {
        HOMER_ANNOTATEPEAKS (
            MACS3_CONSENSUS.out.bed,
            ch_fasta,
            ch_gtf
        )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)

        //
        // MODULE: Add boolean fields to annotated consensus peaks to aid filtering
        //
        ANNOTATE_BOOLEAN_PEAKS (
            MACS3_CONSENSUS.out.boolean_txt.join(HOMER_ANNOTATEPEAKS.out.txt, by: [0])
        )
        ch_versions = ch_versions.mix(ANNOTATE_BOOLEAN_PEAKS.out.versions)
    }

    // Create channels: [ meta, [ ip_bams ], saf ]
    MACS3_CONSENSUS
        .out
        .saf
        .map {
            meta, saf ->
                [ meta.id, meta, saf ]
        }
        .join(ch_bams)
        .map {
            antibody, meta, saf, bams ->
                [ meta, bams.flatten().sort(), saf ]
        }
        .set { ch_bam_saf }

    //
    // Quantify peaks across samples with featureCounts
    //
    SUBREAD_FEATURECOUNTS (
        ch_bam_saf
    )
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect { it -> it[1] })

    //
    // Generate QC plots with DESeq2
    //
    ch_deseq2_qc_pdf           = channel.empty()
    ch_deseq2_qc_rdata         = channel.empty()
    ch_deseq2_qc_rds           = channel.empty()
    ch_deseq2_qc_pca_txt       = channel.empty()
    ch_deseq2_qc_dists_txt     = channel.empty()
    ch_deseq2_qc_log           = channel.empty()
    ch_deseq2_qc_size_factors  = channel.empty()
    if (!skip_deseq2_qc) {
        DESEQ2_QC (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_deseq2_pca_header_multiqc,
            ch_deseq2_clustering_header_multiqc
        )
        ch_deseq2_qc_pdf           = DESEQ2_QC.out.pdf
        ch_deseq2_qc_rdata         = DESEQ2_QC.out.rdata
        ch_deseq2_qc_rds           = DESEQ2_QC.out.rds
        ch_deseq2_qc_pca_txt       = DESEQ2_QC.out.pca_txt
        ch_deseq2_qc_dists_txt     = DESEQ2_QC.out.dists_txt
        ch_deseq2_qc_log           = DESEQ2_QC.out.log
        ch_deseq2_qc_size_factors  = DESEQ2_QC.out.size_factors
        ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.pca_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.dists_multiqc.collect { it -> it[1] })
        ch_versions = ch_versions.mix(DESEQ2_QC.out.versions)
    }

    if (!skip_consensus_plotprofile) {
        MACS3_CONSENSUS
            .out
            .bed
            .map { meta, peaks ->
                [ meta.antibody, meta.exp_type, peaks ]
            }
            .set { ch_cons_peaks }

        ch_bigwigs
            .map { meta, bw ->
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [ antibody, meta.exp_type, meta.norm_factor_type, meta.signal_vs_input_operation, meta.averaged_brep, meta.id, meta, bw ]
            }
            .groupTuple(by: [0, 1, 2, 3, 4])
            // antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, ids, metas, bws
            .combine(ch_cons_peaks, by: [0, 1])
            .map {
                antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, ids, metas, bws, cons_peaks ->
                    def meta_new = metas[0].clone()
                    meta_new.id = exp_type + '_' +
                        (antibody ? antibody : 'no_antibody') +
                        '_' + norm_factor_type +
                        (signal_vs_input_op ? '_' + signal_vs_input_op : '') +
                        (averaged_brep ? '_' + 'bRep_avg' : '')
                    meta_new.antibody = antibody
                    meta_new.ids = ids
                    [ meta_new, bws.flatten(), cons_peaks ]
            }
            .set { ch_bigwigs_peaks }

        //
        // MODULE: deepTools computeMatrix for peaks
        //
        DEEPTOOLS_COMPUTEMATRIX_PEAKS (
            ch_bigwigs_peaks
        )

        //
        // MODULE: deepTools profile plots
        //
        DEEPTOOLS_PLOTPROFILE_PEAKS (
            DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix
        )
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTPROFILE_PEAKS.out.table.collect { it -> it[1] })
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE_PEAKS.out.versions.first())

        //
        // MODULE: deepTools heatmaps
        //
        DEEPTOOLS_PLOTHEATMAP_PEAKS (
            DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP_PEAKS.out.versions.first())
    }


    emit:
    consensus_bed           = MACS3_CONSENSUS.out.bed           // channel: [ bed ]
    consensus_saf           = MACS3_CONSENSUS.out.saf           // channel: [ saf ]
    consensus_pdf           = MACS3_CONSENSUS.out.pdf           // channel: [ pdf ]
    consensus_txt           = MACS3_CONSENSUS.out.txt           // channel: [ pdf ]
    consensus_boolean_txt   = MACS3_CONSENSUS.out.boolean_txt   // channel: [ txt ]
    consensus_intersect_txt = MACS3_CONSENSUS.out.intersect_txt // channel: [ txt ]

    featurecounts_txt       = SUBREAD_FEATURECOUNTS.out.counts  // channel: [ txt ]
    featurecounts_summary   = SUBREAD_FEATURECOUNTS.out.summary // channel: [ txt ]

    deseq2_qc_pdf           = ch_deseq2_qc_pdf                  // channel: [ pdf ]
    deseq2_qc_rdata         = ch_deseq2_qc_rdata                // channel: [ rdata ]
    deseq2_qc_rds           = ch_deseq2_qc_rds                  // channel: [ rds ]
    deseq2_qc_pca_txt       = ch_deseq2_qc_pca_txt              // channel: [ txt ]
    deseq2_qc_dists_txt     = ch_deseq2_qc_dists_txt            // channel: [ txt ]
    deseq2_qc_log           = ch_deseq2_qc_log                  // channel: [ txt ]
    deseq2_qc_size_factors  = ch_deseq2_qc_size_factors         // channel: [ txt ]

    multiqc_files           = ch_multiqc_files                   // channel: [ multiqc_files ]
    versions                = ch_versions                       // channel: [ versions.yml ]
}
