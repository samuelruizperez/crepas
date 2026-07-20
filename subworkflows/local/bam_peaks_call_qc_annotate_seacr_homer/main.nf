
//
// Call peaks with SEACR, annotate with HOMER and perform downstream QC
//

include { SEACR_CALLPEAK } from '../../../modules/nf-core/seacr/callpeak/main'                                                                                                                                                                                                                                                                                  

workflow BAM_PEAKS_CALL_QC_ANNOTATE_SEACR_HOMER {
    
    take:
    ch_bedgraph                            // channel: [ val(meta), [ ip_bedgraph ], [ ipcontrol_bedgraph ] ]
    seacr_peak_threshold

    main:

    SEACR_CALLPEAK (
        ch_bedgraph,
        seacr_peak_threshold
    )
    ch_seacr_peaks       = SEACR_CALLPEAK.out.bed


    emit:
    peaks                        = ch_seacr_peaks                   // channel: [ val(meta), [ peaks ] ]

}
