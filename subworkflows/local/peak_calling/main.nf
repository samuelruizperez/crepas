include { BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER            } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_epic2_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER            } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_macs3_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER          } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_genrich_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER             } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_mace_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER          } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_danpos2_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_consenrich_rocco_homer/main'
include { EDD                                               } from '../../../modules/local/edd/main'
include { DENOPA                                            } from '../../../modules/local/denopa/main'
include { BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2           } from '../../../subworkflows/local/bed_consensus_quantify_qc_bedtools_featurecounts_deseq2/main'

workflow PEAK_CALLING {

    take:
    ch_bam_bai              // channel: [ val(meta), [ bam ], [ bai ] ]
    peak_caller             // String
    ch_fasta
    ch_gtf
    ch_effective_gfraction
    ch_chrom_sizes_endo
    ch_blacklist
    ch_sparsebed
    ch_active_regions
    ch_rocco_params
    ch_effective_gsize
    ch_epic2_peak_count_header
    ch_epic2_frip_score_header
    ch_epic2_peak_annotation_header
    ch_gr_peak_count_header
    ch_gr_frip_score_header
    ch_gr_peak_annotation_header
    ch_mace_peak_count_header
    ch_mace_frip_score_header
    ch_mace_peak_annotation_header
    ch_macs3_peak_count_header
    ch_macs3_frip_score_header
    ch_macs3_peak_annotation_header
    ch_deseq2_pca_header
    ch_deseq2_clustering_header
    narrow_peak
    skip_peak_annotation
    skip_peak_qc
    skip_bdgcmp
    skip_dpeak
    skip_dpos
    skip_consensus_peaks
    skip_deseq2_qc
    skip_consensus_plotprofile
    input_cisrpm_in_plotprofile

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_bam = ch_bam_bai.map { meta, bam, bai -> [meta, bam] }

    //
    // SUBWORKFLOW: Call consensus regions with Consenrich and ROCCO
    //
    ch_consenrich_tracks = channel.empty()
    ch_rocco_peaks = channel.empty()
    if (peak_caller == 'consenrich') {
        BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER (
            ch_bam_bai,
            ch_chrom_sizes_endo,
            ch_blacklist,
            ch_sparsebed.ifEmpty([[:], []]),
            ch_active_regions.ifEmpty([[:], []]),
            ch_rocco_params,
            ch_effective_gsize
        )
        ch_consenrich_tracks = BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.consenrich_signal
        ch_consenrich_tracks = ch_consenrich_tracks.mix(BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.consenrich_residuals)
        ch_consenrich_tracks = ch_consenrich_tracks.mix(BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.consenrich_eratio)
        ch_rocco_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.rocco_peaks
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.versions.first())
    }

    //
    // SUBWORKFLOW: Call peaks with epic2, annotate with HOMER and perform downstream QC
    //
    ch_epic2_peaks = channel.empty()
    ch_epic2_frip_multiqc = channel.empty()
    ch_epic2_peak_count_multiqc = channel.empty()
    ch_epic2_plot_homer_annotatepeaks_tsv = channel.empty()
    if (peak_caller == 'epic2') {
        BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER (
            ch_bam,//.filter { it -> !(it[0].exp_type in ['ChIP-exo', 'OK-seq']) },
            ch_fasta,
            ch_gtf,
            ch_chrom_sizes_endo,
            ch_effective_gfraction,
            ".annotatePeaks.txt",
            ch_epic2_peak_count_header,
            ch_epic2_frip_score_header,
            ch_epic2_peak_annotation_header,
            skip_peak_annotation,
            skip_peak_qc
        )
        ch_epic2_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.peaks
        ch_epic2_frip_multiqc = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.frip_multiqc
        ch_epic2_peak_count_multiqc = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.peak_count_multiqc
        ch_epic2_plot_homer_annotatepeaks_tsv = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.plot_homer_annotatepeaks_tsv
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.versions)
    }

    //
    // SUBWORKFLOW: Call peaks with Genrich, annotate with HOMER and perform downstream QC
    //
    ch_genrich_peaks = channel.empty()
    if (peak_caller == 'genrich') {
        BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER (
            ch_bam,//.filter { it -> !(it[0].exp_type in ['ChIP-exo', 'OK-seq']) },
            ch_fasta,
            ch_gtf,
            ch_blacklist,
            ".annotatePeaks.txt",
            ch_gr_peak_count_header,
            ch_gr_frip_score_header,
            ch_gr_peak_annotation_header,
            narrow_peak,
            skip_peak_annotation,
            skip_peak_qc
        )
        ch_genrich_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.peaks
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.frip_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.peak_count_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it -> it[1] })
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.versions)
    }


    //
    // SUBWORKFLOW: Call peaks with MACE (for ChIP-exo samples)
    //
    ch_mace_peaks = channel.empty()
    if (peak_caller == 'mace') {
        BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER(
            ch_bam_bai,//.filter { it -> it[0].exp_type == 'ChIP-exo' },
            ch_fasta,
            ch_gtf,
            ch_chrom_sizes_endo,
            ".annotatePeaks.txt",
            ch_mace_peak_count_header,
            ch_mace_frip_score_header,
            ch_mace_peak_annotation_header,
            skip_peak_annotation,
            skip_peak_qc
        )
        ch_mace_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.peaks
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.frip_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.peak_count_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it -> it[1] })
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.versions)
    }

    //
    // SUBWORKFLOW: Call peaks with DANPOS2
    //
    if (peak_caller == 'dpeak' || peak_caller == 'dpos') {
        BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER (
            ch_bam,//.filter { it -> !(it[0].exp_type in ['SCAR-seq', 'ChIP-exo', 'OK-seq']) },
            skip_dpeak,
            skip_dpos
        )
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER.out.versions)
    }


    //
    // Create channel for downstream processes: [ meta, [ ip_bam, ipcontrol_bam ] [ ip_bai, ipcontrol_bai ] ]
    // (Excluding ips_wo_ipcontrol as they don't need to be compared to anything)
    //
    ch_bam_bai
        .branch { meta, bam, bai ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, bam, bai]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, bam, bai]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, bam, bai]
        }
        .set { ch_bam_bai_by_type }

    ch_bam_bai_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_bai_by_type.ipcontrols, by: [0,1])
        .map { ipcontrol_id, antibody, ip_meta, ip_bam, ip_bai, ipcontrol_meta, ipcontrol_bam, ipcontrol_bai ->
            [ ip_meta, [ip_bam] + [ipcontrol_bam], [ip_bai] + [ipcontrol_bai] ]
        }
        .set { ch_ip_and_ipcontrols_bam_bai }

    // Create channels: [ meta, ip_bam, ipcontrol_bam ]
    // Including ips_wo_ipcontrol as they will be used for peak calling without control
    ch_bam_bai_by_type
        .ips_wo_ipcontrol
        .map { meta, bam, bai -> [meta, [bam], [bai]] }
        .mix(ch_ip_and_ipcontrols_bam_bai)
        // ips_wo_ipcontrol do not have ipcontrol_bam
        .map { meta, bams, bais ->
            [meta, bams[0], (bams[1] ?: [])]
        }
        .set { ch_all_ip_and_controls }


    ch_edd_peaks = channel.empty()
    if (peak_caller == 'edd') {
        //
        // MODULE: Call peaks with EDD
        //
        EDD (
            ch_all_ip_and_controls,
            ch_chrom_sizes_endo,
            ch_blacklist
        )
        ch_edd_peaks = EDD.out.peaks
        ch_versions = ch_versions.mix(EDD.out.versions.first())
    }

    ch_denopa_peaks = channel.empty()
    if (peak_caller == 'denopa') {
        //
        // MODULE: Call peaks with denopa
        //
        DENOPA (
            ch_all_ip_and_controls.filter { it -> it[0].exp_type in ['ATAC-seq'] }
        )
        ch_denopa_peaks = DENOPA.out.arers
        ch_versions = ch_versions.mix(DENOPA.out.versions.first())
    }

    ch_macs3_peaks = channel.empty()
    if (peak_caller == 'macs3') {
        //
        // SUBWORKFLOW: Call peaks with MACS3, annotate with HOMER and perform downstream QC
        //
        BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER (
            ch_all_ip_and_controls,
            ch_fasta,
            ch_gtf,
            ch_chrom_sizes_endo,
            ch_effective_gsize,
            "_peaks.annotatePeaks.txt",
            ch_macs3_peak_count_header,
            ch_macs3_frip_score_header,
            ch_macs3_peak_annotation_header,
            narrow_peak,
            skip_peak_annotation,
            skip_peak_qc,
            skip_bdgcmp
        )
        ch_macs3_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peaks
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.frip_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peak_count_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it -> it[1] })
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.versions)
    }


    //
    //  Consensus peaks analysis
    //
    ch_consensus_bed = channel.empty()
    ch_consensus_txt = channel.empty()
    if (!skip_consensus_peaks) {
        // Create channels: [ antibody, [ ip_bams ] ]
        ch_all_ip_and_controls
            .map { meta, ip_bam, ipcontrol_bam ->
                [meta.antibody, ip_bam]
            }
            .groupTuple()
            .set { ch_antibody_bams }

        BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2(
            ch_macs3_peaks,
            ch_antibody_bams,
            BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_all_endo,
            ch_fasta.map { it -> it[1] },
            ch_gtf.map { it -> it[1] },
            ch_deseq2_pca_header,
            ch_deseq2_clustering_header,
            narrow_peak,
            skip_peak_annotation,
            skip_deseq2_qc,
            skip_consensus_plotprofile,
            input_cisrpm_in_plotprofile
        )
        ch_consensus_bed = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.consensus_bed
        ch_consensus_txt = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.consensus_txt
        ch_multiqc_files = ch_multiqc_files.mix(BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.multiqc_files)
        ch_versions = ch_versions.mix(BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.versions)
    }

    emit:

    //bedgraph_endo    = ch_bdg_all.filter { it -> it[0].genome == genome }      // channel: [ val(meta), [ bedgraph ] ]

    versions      = ch_versions                                     // channel: [ versions.yml ]
}
