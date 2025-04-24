include { DEEPTOOLS_BAMCOVERAGE }       from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { BEDTOOLS_MAKEWINDOWS }        from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAP }                from '../../../modules/nf-core/bedtools/map/main'
include { BEDGRAPH_NORMALIZE }          from '../../../modules/local/bedgraph_normalize/main'
include { BEDGRAPH_SIGNAL_OVER_INPUT }  from '../../../modules/local/bedgraph_signal_over_input/main'
include { UCSC_BEDGRAPHTOBIGWIG     }   from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BAM_NORMALIZED_BIGWIG_DEEPTOOLS {

    take:
    ch_bam_bai               // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_chrom_sizes          // channel: [ bed ]
    genome                  // string: genome name
    spikein_genome          // string: spike-in genome name
    skip_srpm           // boolean: skip the SRPM normalization step
    skip_cisrpm          // boolean: skip the CISRPM normalization step
    skip_cisrpmsoi      // boolean: skip the CISRPM-SOI normalization step

    main:

    ch_versions = Channel.empty()


    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1
            meta_clone.norm_factor_type = 'raw'
            [ meta_clone, bam, bai ]
        }
        .set { ch_bam_bai_raw }

    //
    // MODULE: Calculate raw coverage per bin
    //
    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai_raw,
        [],
        []
    )
    ch_bdg_raw = DEEPTOOLS_BAMCOVERAGE.out.bedgraph
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())

    // bamCoverage merges contiguous bins with the same coverage and there is no option to disable this.
    // See https://github.com/deeptools/deepTools/issues/907#issuecomment-576729674
    // Therefore, to get bins of equal size across the whole genome,
    // we have to make windows and map the bedgraph to the windows:

    //
    // MODULE: Make windows of equal size across the whole genome
    //
    BEDTOOLS_MAKEWINDOWS (
        ch_chrom_sizes
    )
    ch_windows = BEDTOOLS_MAKEWINDOWS.out.bed
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())

    // Create channel: [ val(meta_bdg_raw), windows, bdg_raw ]
    ch_bdg_raw
        .combine(ch_windows)
        .map { meta_bdg_raw, bdg_raw, meta_windows, windows ->
            [ meta_bdg_raw, windows, bdg_raw ]
        }
        .set { ch_windows_bdg_raw }

    //
    // MODULE: Map the coverage bedgraph to the windows
    //
    BEDTOOLS_MAP (
        ch_windows_bdg_raw,
        ch_chrom_sizes      
    )
    ch_bdg_raw = BEDTOOLS_MAP.out.bed
    ch_versions = ch_versions.mix(BEDTOOLS_MAP.out.versions.first())

    // Modify channel meta to add RPM normalization factors
    ch_bdg_rpm = Channel.empty()

    ch_bdg_raw
        .map {
            meta, bdg ->
                def meta_clone = meta.clone()
                meta_clone.norm_factor_val = 1e6 / meta_clone.total_mapped_reads
                meta_clone.norm_factor_type = 'rpm'
                [ meta_clone, bdg ]
        }
        .set { ch_bdg_rpm }

    // Modify channel meta to add RPM normalization factors
    ch_bdg_srpm = Channel.empty()
    if (!skip_srpm) {

        ch_bdg_raw
            .map {
                meta, bdg ->
                    [ meta.id, meta, bdg ]
            }
            .branch { id, meta, bdg ->
                endo: meta.genome == genome
                exo: meta.genome == spikein_genome
            }
            .set { ch_bdg_genome }

        ch_bdg_genome
            .endo
            .combine(ch_bdg_genome.exo, by: 0)
            .map { id, endo_meta, endo_bdg, exo_meta, exo_bdg ->
                def meta_clone = endo_meta.clone()
                meta_clone.norm_factor_val = 1e6 / exo_meta.total_mapped_reads
                meta_clone.norm_factor_type = 'srpm'
                [ meta_clone, endo_bdg ]
            }
            .set { ch_bdg_srpm }
    }

    ch_bdg_cisrpm = Channel.empty()
    if (!skip_cisrpm) {

        // Split BAMs by genome (endo and exo) and by type (ip and control)
        ch_bdg_raw
            .map {
                meta, bdg ->
                    [ meta.id, meta, bdg ]
            }
            .branch { id, meta, bdg ->
                endo_ip: meta.genome == genome && !meta.is_control
                endo_control: meta.genome == genome && meta.is_control
                exo_ip: meta.genome == spikein_genome && !meta.is_control
                exo_control: meta.genome == spikein_genome && meta.is_control
            }
            .set { ch_bdg_genome_type }

        // Combine the endo and exo BAMs (ChIPs)
        ch_bdg_genome_type
            .endo_ip
            .combine(ch_bdg_genome_type.exo_ip, by: 0)
            .map { 
                ip_id, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ->
                    [ endo_ip_meta.control, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ]
            }
            .set { ch_bdg_genome_ip }

        // Combine the endo and exo BAMs (inputs)
        ch_bdg_genome_type
            .endo_control
            .combine(ch_bdg_genome_type.exo_control, by: 0)
            .map { control_id, endo_control_meta, endo_control_bdg, exo_control_meta, exo_control_bdg ->
                [ endo_control_meta.id, endo_control_meta, endo_control_bdg, exo_control_meta, exo_control_bdg ]
            }
            .set { ch_bdg_genome_control }

        // Combine the combined ChIPs with the combined inputs
        ch_bdg_genome_ip
            .combine(ch_bdg_genome_control, by: 0)
            .map { 
                id, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg, endo_control_meta, endo_control_bdg, exo_control_meta, exo_control_bdg ->
                    def meta_clone = endo_ip_meta.clone()
                    meta_clone.norm_factor_val = (1e6 / exo_ip_meta.total_mapped_reads) * (exo_control_meta.total_mapped_reads / endo_control_meta.total_mapped_reads)
                    meta_clone.norm_factor_type = 'cisrpm'
                    [ meta_clone, endo_ip_bdg ]
            }
            .set { ch_bdg_ip_cisrpm }
        
        // Now do the missing CISRPM for the endogenous inputs
        // In this case CISRPM is the same as RPM
        ch_bdg_rpm
            .filter { meta, bdg ->
                meta.genome == genome && meta.is_control
            }
            .map { meta, bdg ->
                def meta_clone = meta.clone()
                // just change norm_factor_type to cisrpm but keeping norm_factor_val the untouched
                meta_clone.norm_factor_type = 'cisrpm'
                [ meta_clone, bdg ]
            }
            .set { ch_bdg_control_cisrpm }

        ch_bdg_cisrpm = ch_bdg_ip_cisrpm.mix(ch_bdg_control_cisrpm)
    }

    ch_bdg_rpm
        .mix(ch_bdg_srpm)
        .mix(ch_bdg_cisrpm)
        .set { ch_bdg_norm }

    // TODO: print for debugging
    ch_bdg_norm
        .map {
            meta, bdg ->
                "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_norm.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    //
    // MODULE: Normalize the bedgraphs
    //
    BEDGRAPH_NORMALIZE (
        ch_bdg_norm
    )
    ch_versions = ch_versions.mix(BEDGRAPH_NORMALIZE.out.versions.first())


    // Add raw bedgraph to 
    ch_bdg_all = ch_bdg_raw.mix(BEDGRAPH_NORMALIZE.out.bedgraph)

    // //
    // // MODULE: Sort the bedgraph
    // //
    // FILE_SORT_NORM (
    //    DEEPTOOLS_BAMCOVERAGE.out.bedgraph,
    //     'bedgraph'
    // )
    // ch_bdg = FILE_SORT_NORM.out.sorted
    // ch_versions = ch_versions.mix(FILE_SORT_NORM.out.versions.first())


    // Create channel: [ val(meta), [ ip_bdg ], [ control_bdg ] ]
    ch_bdg_ip_control_cisrpm = Channel.empty()
    if (!skip_cisrpmsoi) {
        ch_bdg_all
            .filter { meta, bdg ->
                meta.norm_factor_type == 'cisrpm'
            }
            .branch { meta, bdg ->
                ips_with_control: meta.control
                    return [ meta.control, meta, bdg ]
                // Cannot calculate CISRPM-SOI for ChIPs without inputs
                // ips_without_control: !meta.control && !meta.is_control
                //     return [ meta, bdg ]
                controls: !meta.control && meta.is_control
                    return [ meta.id, meta, bdg ]
            }
            .set { ch_bdg_ip_control_cisrpm }

        ch_bdg_ip_control_cisrpm
            .ips_with_control
            .combine(ch_bdg_ip_control_cisrpm.controls, by: 0)
            .map { control_id, ip_meta, ip_bdg, control_meta, control_bdg ->
                def meta_clone = ip_meta.clone()
                meta_clone.signal_over_input = true
                [ meta_clone, ip_bdg, control_bdg ]
            }
            .set { ch_bdg_ip_control_cisrpm }

        //
        // MODULE: Calculate CIRSPM signal over input (ChIP over input)
        //
        BEDGRAPH_SIGNAL_OVER_INPUT (
            ch_bdg_ip_control_cisrpm
        )
        ch_versions = ch_versions.mix(BEDGRAPH_SIGNAL_OVER_INPUT.out.versions.first())

        ch_bdg_all = BEDGRAPH_SIGNAL_OVER_INPUT.out.bedgraph.mix(ch_bdg_all)

        // //
        // // MODULE: Sort the SOI bedgraph
        // //
        // FILE_SORT_SOI (
        //     BEDGRAPH_SIGNAL_OVER_INPUT.out.bedgraph,
        //     'bedgraph'
        // )

        // ch_bdg_all = FILE_SORT_SOI.out.sorted.mix(ch_bdg_all)

    }

    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG (
        ch_bdg_all,
        ch_chrom_sizes.map { it[1] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())
   

    emit:
    bigwig      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig        // channel: [ val(meta), [ bigwig ] ]

    versions      = ch_versions                            // channel: [ versions.yml ]
}

