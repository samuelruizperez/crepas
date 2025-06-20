include { FILE_SORT as SIZES_SORT } from '../../../modules/local/file_sort/main'
include { CONSENRICH           } from '../../../modules/local/consenrich/main'
include { ROCCO               } from '../../../modules/local/rocco/main'


workflow BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_HOMER {
    take:
    ch_bam                // channel: [ meta, [ip_bams_merged_reps], [ip_bais_merged_reps], [control_bams_merged_reps], [control_bais_merged_reps] ]
    ch_chrom_sizes        // channel: [ val(meta), [ chrom_sizes ] ]
    ch_blacklist          // channel: [ val(meta), [ blacklist ] ]
    ch_sparsebed          // channel: [ val(meta), [ sparsebed ] ]
    ch_active_regions     // channel: [ val(meta), [ active_regions ] ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Sort chromosome sizes (avoid Consenrich errors)
    //
    SIZES_SORT (
        ch_chrom_sizes,
        'sizes'
    )
    ch_chrom_sizes = SIZES_SORT.out.sorted
    ch_versions = ch_versions.mix(SIZES_SORT.out.versions.first())


    //
    // MODULE: Integrate ChIPs and their input controls
    //
    CONSENRICH (
        ch_bam,
        ch_chrom_sizes,
        ch_blacklist,
        ch_sparsebed,
        ch_active_regions
    )
    ch_versions = ch_versions.mix(CONSENRICH.out.versions.first())

    
    // Create channel: [ val(meta), ch_csr_signal, bamlist_txt ]
    CONSENRICH
        .out
        .signal_track
        .combine(ch_bam, by: 0)
        // this is: [ val(meta), ch_csr_signal, [ ip_bams ], [ ip_bais ], [ control_bams ], [ control_bais ] ]
        .map { meta, ch_csr_signal, ip_bams, ip_bais, control_bams, control_bais ->
            [ ]




    //
    // MODULE: Call peaks
    //
    ROCCO (
        CONSENRICH.out.results,
        ch_bam_list,
        ch_chrom_sizes,
        CONSENRICH.out.params_file,
        CONSENRICH.out.effective_genome_size
    )
    ch_versions = ch_versions.mix(ROCCO.out.versions.first())

    emit:

    results     = CONSENRICH.out.results   // channel: [ val(meta), [ tsv ] ]

    versions    = ch_versions              // channel: [ versions.yml ]
}
