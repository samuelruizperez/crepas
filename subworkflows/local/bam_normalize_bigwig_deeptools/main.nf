include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_BINS           } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_ENDO             } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_EXO              } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_ENDO                             } from '../../../modules/nf-core/bedtools/map/main'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_EXO                              } from '../../../modules/nf-core/bedtools/map/main'
include { BEDGRAPH_NORMALIZE                                            } from '../../../modules/local/bedgraph_normalize/main'
include { BEDGRAPH_SIGNAL_OVER_INPUT                                    } from '../../../modules/local/bedgraph_signal_over_input/main'
include { FILE_SORT as BEDGRAPH_SORT                                    } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_ENDO   } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_EXO    } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { DEEPTOOLS_BIGWIGCOMPARE } from '../../../modules/nf-core/deeptools/bigwigcompare/main'
include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_BINSIZE1       } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { DEEPTOOLS_BIGWIGAVERAGE as DEEPTOOLS_BIGWIGAVERAGE_BINS        } from '../../../modules/local/deeptools/bigwigaverage/main'
include { DEEPTOOLS_BIGWIGAVERAGE as DEEPTOOLS_BIGWIGAVERAGE_BINSIZE1     } from '../../../modules/local/deeptools/bigwigaverage/main'

workflow BAM_NORMALIZE_BIGWIG_DEEPTOOLS {

    take:
    ch_bam_bai              // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_chrom_sizes_endo     // channel: [ bed ]
    ch_chrom_sizes_exo      // channel: [ bed ]
    coverage_bin_size       // int: size of the coverage bin in bp
    genome                  // string: genome name
    spikein_genome          // string: spike-in genome name
    skip_srpm               // boolean: skip the SRPM normalization step
    skip_cisrpm             // boolean: skip the CISRPM normalization step
    skip_cisrpmsoi          // boolean: skip the CISRPM-SOI normalization step
    skip_plot_profile       // boolean: skip the plot profile step
    min_reads_for_norm
    skip_rpm_lfc
    skip_bw_average

    main:

    ch_versions = channel.empty()


    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1
            meta_clone.norm_factor_type = 'raw'
            [ meta_clone, bam, bai ]
        }
        // Remove empty BAMs to prevent bamCoverage errors
        .filter { meta, bam, bai -> meta[meta.last_total_mapped_reads_key] >= min_reads_for_norm }
        .set { ch_bam_bai }

    //
    // MODULE: Calculate raw coverage per bin
    //
    DEEPTOOLS_BAMCOVERAGE_BINS (
        ch_bam_bai,
        [],
        []
    )
    ch_bdg_raw = DEEPTOOLS_BAMCOVERAGE_BINS.out.bedgraph
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_BINS.out.versions.first())

    // TODO: print for debugging
    ch_bdg_raw
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_raw.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    // bamCoverage merges contiguous bins with the same coverage and there is no option to disable this.
    // See https://github.com/deeptools/deepTools/issues/907#issuecomment-576729674
    // Therefore, to get bins of equal size across the whole genome,
    // we have to make windows and map the bedgraph to the windows.
    // An alternative would be to use featureCounts:

    //
    // MODULE: Make windows of equal size across the whole genome
    //
    BEDTOOLS_MAKEWINDOWS_ENDO (
        ch_chrom_sizes_endo
    )
    ch_windows_endo = BEDTOOLS_MAKEWINDOWS_ENDO.out.bed
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS_ENDO.out.versions)

    // Create channel: [ val(meta_bdg_raw), windows, bdg_raw ]
    ch_bdg_raw
        .filter { meta, bdg ->
            meta.genome == genome
        }
        .combine(ch_windows_endo)
        .map { meta_bdg_raw, bdg_raw, meta_windows, windows ->
            [ meta_bdg_raw, windows, bdg_raw ]
        }
        .set { ch_windows_endo_bdg_raw }

    //
    // MODULE: Map the coverage bedgraph to the windows (endogenous genome)
    //
    BEDTOOLS_MAP_ENDO (
        ch_windows_endo_bdg_raw,
        ch_chrom_sizes_endo
    )
    ch_bdg_map_endo = BEDTOOLS_MAP_ENDO.out.mapped
    ch_versions = ch_versions.mix(BEDTOOLS_MAP_ENDO.out.versions.first())
    
    ch_bdg_map = ch_bdg_map_endo
    ch_windows_exo = channel.empty()
    ch_windows_exo_bdg_raw = channel.empty()
    if (spikein_genome) {
        BEDTOOLS_MAKEWINDOWS_EXO (
            ch_chrom_sizes_exo
        )
        ch_windows_exo = BEDTOOLS_MAKEWINDOWS_EXO.out.bed
        ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS_EXO.out.versions)

        // Create channel: [ val(meta_bdg_raw), windows, bdg_raw ]
        ch_bdg_raw
            .filter { meta, bdg ->
                meta.genome == spikein_genome
            }
            .combine(ch_windows_exo)
            .map { meta_bdg_raw, bdg_raw, meta_windows, windows ->
                [ meta_bdg_raw, windows, bdg_raw ]
            }
            .set { ch_windows_exo_bdg_raw }

        //
        // MODULE: Map the coverage bedgraph to the windows (spike-in genome)
        //
        BEDTOOLS_MAP_EXO (
            ch_windows_exo_bdg_raw,
            ch_chrom_sizes_exo
        )
        ch_bdg_map_exo = BEDTOOLS_MAP_EXO.out.mapped
        ch_versions = ch_versions.mix(BEDTOOLS_MAP_EXO.out.versions.first())

        // Merge the two channels
        ch_bdg_map = ch_bdg_map_endo.mix(ch_bdg_map_exo)
    }

    // TODO: print for debugging
    ch_bdg_map
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_map.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    // RPM normalization factors
    ch_bdg_rpm = channel.empty()
    ch_bdg_map
        .map { meta, bdg ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1e6 / meta[meta.ref_total_mapped_reads_for_rpm]
            meta_clone.norm_factor_type = 'rpm'
            [ meta_clone, bdg ]
        }
        .set { ch_bdg_rpm }

    // TODO: print for debugging
    ch_bdg_rpm
        .map {
            meta, bdg ->
                "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_rpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    // Copy and modify channel meta to add SRPM normalization factors (for ChIPs)
    ch_bdg_srpm = channel.empty()
    if (!skip_srpm) {
        ch_bdg_map
            .map { meta, bdg ->
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody_to_use = meta.antibody ?: meta.input_control_of_antibody
                [ meta.id, antibody_to_use, meta, bdg ]
            }
            .branch { id, antibody, meta, bdg ->
                endo: meta.genome == genome
                exo: meta.genome == spikein_genome
            }
            .set { ch_bdg_genome }

        // TODO: print for debugging
        ch_bdg_genome.endo
            .map { it ->
                "${it}"
            }
            .collectFile( name: 'ch_bdg_srpm_genome_endo.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

        // TODO: print for debugging
        ch_bdg_genome.exo
            .map { it ->
                "${it}"
            }
            .collectFile( name: 'ch_bdg_srpm_genome_exo.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )


        ch_bdg_genome.endo
            .combine(ch_bdg_genome.exo, by: [0,1])
            .map { id, antibody, endo_meta, endo_bdg, exo_meta, exo_bdg ->
                def meta_clone = endo_meta.clone()
                meta_clone.norm_factor_val = 1e6 / exo_meta[exo_meta.ref_total_mapped_reads_for_srpm]
                meta_clone.norm_factor_type = 'srpm'
                [ meta_clone, endo_bdg ]
            }
            .set { ch_bdg_srpm }

    }


    // TODO: print for debugging
    ch_bdg_srpm
        .map { it ->
            "${it}"
        }
        .collectFile( name: 'ch_bdg_srpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    
    // Copy and modify channel meta to add CISRPM normalization factors
    ch_bdg_genome_type = channel.empty()
    ch_bdg_genome_ip = channel.empty()
    ch_bdg_genome_ipcontrol = channel.empty()
    ch_bdg_ip_cisrpm = channel.empty()
    ch_bdg_ipcontrol_cisrpm = channel.empty()
    ch_bdg_cisrpm = channel.empty()
    // "if (spikein_genome)" is needed, otherwise cisrpm will be attempted for
    // controls, and this will fail, since there is no ref_total_mapped_reads_for_cisrpm
    if (spikein_genome && !skip_cisrpm) {
        // Split BAMs by genome (endo and exo) and by type (ip and ipcontrol)
        ch_bdg_map
            .map { meta, bdg ->
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody_to_use = meta.antibody ?: meta.input_control_of_antibody
                [ meta.id, antibody_to_use, meta, bdg ]
            }
            .branch { id, antibody, meta, bdg ->
                endo_ip: meta.genome == genome && !meta.is_input_control
                endo_ipcontrol: meta.genome == genome && meta.is_input_control
                exo_ip: meta.genome == spikein_genome && !meta.is_input_control
                exo_ipcontrol: meta.genome == spikein_genome && meta.is_input_control
            }
            .set { ch_bdg_genome_type }

        // Combine the endo and exo BAMs (ChIPs)
        ch_bdg_genome_type.endo_ip
            .combine(ch_bdg_genome_type.exo_ip, by: [0,1])
            .map { ip_id, ip_antibody, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ->
                    [ endo_ip_meta.input_control, endo_ip_meta.antibody, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ]
            }
            .set { ch_bdg_genome_ip }

        // TODO: print for debugging
        ch_bdg_genome_ip
            .map {
                ipcontrol_id, ip_antibody, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ->
                    "${ipcontrol_id}\t${ip_antibody}\t${endo_ip_meta}\t${endo_ip_bdg}\t${exo_ip_meta}\t${exo_ip_bdg}"
            }
            .collectFile( name: 'ch_bdg_genome_ip.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

        // Combine the endo and exo BAMs (inputs)
        ch_bdg_genome_type.endo_ipcontrol
            .combine(ch_bdg_genome_type.exo_ipcontrol, by: [0,1])
            .map { ipcontrol_id, ipcontrol_antibody, endo_ipcontrol_meta, endo_ipcontrol_bdg, exo_ipcontrol_meta, exo_ipcontrol_bdg ->
                [ endo_ipcontrol_meta.id, endo_ipcontrol_meta.input_control_of_antibody, endo_ipcontrol_meta, endo_ipcontrol_bdg, exo_ipcontrol_meta, exo_ipcontrol_bdg ]
            }
            .set { ch_bdg_genome_ipcontrol }
        
        // TODO: print for debugging
        ch_bdg_genome_ipcontrol
            .map {
                ipcontrol_id, ipcontrol_antibody, endo_ipcontrol_meta, endo_ipcontrol_bdg, exo_ipcontrol_meta, exo_ipcontrol_bdg ->
                    "${ipcontrol_id}\t${ipcontrol_antibody}\t${endo_ipcontrol_meta}\t${endo_ipcontrol_bdg}\t${exo_ipcontrol_meta}\t${exo_ipcontrol_bdg}"
            }
            .collectFile( name: 'ch_bdg_genome_ipcontrol.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

        // Combine the combined ChIPs with the combined inputs
        ch_bdg_genome_ip
            .combine(ch_bdg_genome_ipcontrol, by: [0,1])
            .map { id, antibody, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg, endo_ipcontrol_meta, endo_ipcontrol_bdg, exo_ipcontrol_meta, exo_ipcontrol_bdg ->
                    def meta_clone = endo_ip_meta.clone()
                    meta_clone.norm_factor_val = (1e6 / exo_ip_meta[exo_ip_meta.ref_total_mapped_reads_for_cisrpm]) * (exo_ipcontrol_meta[exo_ipcontrol_meta.ref_total_mapped_reads_for_cisrpm] / endo_ipcontrol_meta[endo_ipcontrol_meta.ref_total_mapped_reads_for_cisrpm])
                    meta_clone.norm_factor_type = 'cisrpm'
                    [ meta_clone, endo_ip_bdg ]
            }
            .set { ch_bdg_ip_cisrpm }

        // TODO: print for debugging
        ch_bdg_ip_cisrpm
            .map {
                meta, bdg ->
                    "${meta}\t${bdg}"
            }
            .collectFile( name: 'ch_bdg_ip_cisrpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )
        
        // Now do the missing CISRPM for the endogenous inputs
        // In this case CISRPM is the same as RPM, but we cannot
        // just copy the RPM from before, since ref_total_mapped_reads_for_cisrpm
        // can be different than ref_total_mapped_reads_for_rpm
        ch_bdg_map
            .filter { meta, bdg ->
                meta.genome == genome && meta.is_input_control
            }
            .map { meta, bdg ->
                def meta_clone = meta.clone()
                meta_clone.norm_factor_val = 1e6 / meta[meta.ref_total_mapped_reads_for_cisrpm]
                meta_clone.norm_factor_type = 'cisrpm'
                [ meta_clone, bdg ]
            }
            .set { ch_bdg_ipcontrol_cisrpm }


        // TODO: print for debugging
        ch_bdg_ipcontrol_cisrpm
            .map {
                meta, bdg ->
                    "${meta}\t${bdg}"
            }
            .collectFile( name: 'ch_bdg_ipcontrol_cisrpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

        ch_bdg_cisrpm = ch_bdg_ip_cisrpm.mix(ch_bdg_ipcontrol_cisrpm)
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
        .collectFile( name: 'ch_bdg_norm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )


    //
    // MODULE: Normalize the bedgraphs
    //
    BEDGRAPH_NORMALIZE (
        ch_bdg_norm,
        'bedGraph'
    )
    ch_versions = ch_versions.mix(BEDGRAPH_NORMALIZE.out.versions.first())

    ch_bdg_map_norm = ch_bdg_map.mix(BEDGRAPH_NORMALIZE.out.normalized)

    // Create channel: [ val(meta), [ ip_bdg ], [ ipcontrol_bdg ] ]
    ch_bdg_all = ch_bdg_map_norm
    ch_bdg_ip_control_cisrpm = channel.empty()
    if (!skip_cisrpmsoi) {
        ch_bdg_map_norm
            .filter { meta, bdg ->
                meta.norm_factor_type == 'cisrpm'
            }
            .branch { meta, bdg ->
                ips_with_ipcontrol: meta.input_control
                    return [ meta.input_control, meta.antibody, meta, bdg ]
                // Cannot calculate CISRPM-SOI for ChIPs without inputs
                // ips_without_ipcontrol: !meta.input_control && !meta.is_input_control
                //     return [ meta, bdg ]
                ipcontrols: !meta.input_control && meta.is_input_control
                    return [ meta.id, meta.input_control_of_antibody, bdg ]
            }
            .set { ch_bdg_ip_control_cisrpm }

        ch_bdg_ip_control_cisrpm
            .ips_with_ipcontrol
            .combine(ch_bdg_ip_control_cisrpm.ipcontrols, by: [0,1])
            .map { ipcontrol_id, ip_antibody, ip_meta, ip_bdg, ipcontrol_bdg ->
                def meta_clone = ip_meta.clone()
                    meta_clone.signal_over_input = true
                    [ meta_clone, ip_bdg, ipcontrol_bdg ]
            }
            .set { ch_bdg_ip_control_cisrpm }


        //
        // MODULE: Calculate CIRSPM signal over input (ChIP over input)
        //
        BEDGRAPH_SIGNAL_OVER_INPUT (
            ch_bdg_ip_control_cisrpm
        )
        ch_bdg_all = BEDGRAPH_SIGNAL_OVER_INPUT.out.bedgraph.mix(ch_bdg_map_norm)
        ch_versions = ch_versions.mix(BEDGRAPH_SIGNAL_OVER_INPUT.out.versions.first())

    }


    //
    // MODULE: Sort the bedgraph so that it works with ucsc_bedgraphtobigwig
    //
    BEDGRAPH_SORT (
        ch_bdg_all,
        'bedGraph'
    )
    ch_bdg_all = BEDGRAPH_SORT.out.sorted
    ch_versions = ch_versions.mix(BEDGRAPH_SORT.out.versions.first())

  
    // Remove empty bedgraphs (mostly from SOI) to avoid bdg to bw conversion errors
    ch_bdg_all
        .filter {
            meta, bdg ->
                bdg.size() > 0
        }
        .set { ch_bdg_all }

    // TODO: print for debugging
    ch_bdg_all
        .map {
            meta, bdg ->
                "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_all.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )


    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_ENDO (
        ch_bdg_all.filter { it -> it[0].genome == genome },
        ch_chrom_sizes_endo.map { it -> it[1] }
    )
    ch_bigwig_endo_rpm = UCSC_BEDGRAPHTOBIGWIG_ENDO.out.bigwig.filter { it -> it[0].norm_factor_type == 'rpm' }
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_ENDO.out.versions.first())

    //
    // MODULE: Compute log2FC between engogenous ChIPs and endogenous inputs
    //
    ch_bigwig_endo_rpm_lfc = channel.empty()
    if (!skip_rpm_lfc) {

        // Create channel: [ val(meta), [ ip_bigwig ], [ ipcontrol_bigwig ] ]
        ch_bigwig_endo_rpm
            .branch { meta, bw ->
                ips: !meta.is_input_control
                    return [ meta.input_control, meta.antibody, meta, bw ]
                ipcontrols: meta.is_input_control
                    return [ meta.id, meta.input_control_of_antibody, meta, bw ]
            }
            .set { ch_bigwig_endo_rpm_type }

        ch_bigwig_endo_rpm_type
            .ips
            .combine(ch_bigwig_endo_rpm_type.ipcontrols, by: [0,1])
            .map { id, antibody, ip_meta, ip_bw, ipcontrol_meta, ipcontrol_bw ->
                def meta_clone = ip_meta.clone()
                    meta_clone.log2FC = true
                    [ meta_clone, ip_bw, ipcontrol_bw ]
            }
            .set { ch_bigwig_endo_rpm_ip_ipcontrol }

        DEEPTOOLS_BIGWIGCOMPARE (
            ch_bigwig_endo_rpm_ip_ipcontrol,
            channel.value([[:], []])
        )
        ch_bigwig_endo_rpm_lfc = DEEPTOOLS_BIGWIGCOMPARE.out.output
        ch_versions = ch_versions.mix(DEEPTOOLS_BIGWIGCOMPARE.out.versions.first())
    }

    ch_bw_exo = channel.empty()
    if (spikein_genome) {
        UCSC_BEDGRAPHTOBIGWIG_EXO (
            ch_bdg_all.filter { it -> it[0].genome == spikein_genome },
            ch_chrom_sizes_exo.map { it -> it[1] }
        )
        ch_bw_exo = UCSC_BEDGRAPHTOBIGWIG_EXO.out.bigwig
        ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_EXO.out.versions.first())
    }

    ch_bw_avg_endo = channel.empty()
    if (!skip_bw_average) {

        // Create channel: [ val(meta), [ bRep_bigwigs ] ]
        UCSC_BEDGRAPHTOBIGWIG_ENDO
            .out
            .bigwig
            .map { meta, bw ->
                def meta_clone = meta.clone()
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                meta_clone.id = meta_clone.id - ~/_bRep_.*$/
                [ meta_clone.id, antibody, meta_clone, bw ]
            }
            .groupTuple(by: [0, 1])
            .map { id, antibody, metas, bws ->
                def meta_clone = metas[0].clone()
                meta_clone.averaged_brep = true
                [meta_clone, bws.flatten()]
            }
            .set { ch_bdg_all_brep_bw }

        //
        // MODULE: Average bigwigs across replicates
        //
        DEEPTOOLS_BIGWIGAVERAGE_BINS (
            ch_bdg_all_brep_bw,
            channel.value([[:], []])
        )
        ch_bw_avg_endo = DEEPTOOLS_BIGWIGAVERAGE_BINS.out.bigwig
        ch_versions = ch_versions.mix(DEEPTOOLS_BIGWIGAVERAGE_BINS.out.versions.first())
    }

    // if coverage_bin_size is not 1, then we need to generate bw with that binsize for computeMatrix
    ch_binsize1 = channel.empty()
    ch_bw_avg_binsize1 = channel.empty()
    if (coverage_bin_size != 1 && !skip_plot_profile) {

        ch_bam_bai
            .branch { meta, bam, bai ->
                ip: !meta.is_input_control
                ipcontrol: meta.is_input_control
            }
            .set { ch_bam_bai_type }

        ch_bam_bai_type
            // First, modify the controls' metas to add their corresponding ChIP's antibody
            .ipcontrol
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .combine(ch_bam_bai_type.ip.map { meta, bam, bai -> [ meta.input_control, meta, bam, bai ] }, by: 0)
            // Temporarily put the meta.antibody in the ipcontrol meta
            .map { ipcontrol_id, ipcontrol_meta, ipcontrol_bam, ipcontrol_bai, ip_meta, ip_bam, ip_bai ->
                    def meta_clone = ipcontrol_meta.clone()
                    meta_clone.input_control_of_antibody = ip_meta.antibody
                    [ meta_clone, ipcontrol_bam, ipcontrol_bai ]
            }
            // remove duplicates based on ipcontrol_meta and filename (basically the ipcontrol_meta.antibody we added above)
            // Because we don't need the same input normalized in the same way several times
            .unique()
            // Now we mix the ipcontrol and ip channels to evaluate them together below
            .mix(ch_bam_bai_type.ip)
            .set { ch_bam_bai_mod }

        ch_bam_bai_mod
            // we only want these bw for the endo genome
            .filter { meta, bam, bai ->
                meta.genome == genome
            }
            .map { meta, bam, bai ->
                def meta_clone = meta.clone()
                meta_clone.norm_factor_val = 1e6 / meta[meta.ref_total_mapped_reads_for_rpm]
                meta_clone.norm_factor_type = 'rpm'
                    [ meta_clone, bam, bai ]
            }
            .set { ch_binsize1 }
        
    
        DEEPTOOLS_BAMCOVERAGE_BINSIZE1 (
            ch_binsize1,
            [],
            []
        )
        ch_binsize1 = DEEPTOOLS_BAMCOVERAGE_BINSIZE1.out.bigwig
        ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_BINSIZE1.out.versions.first())

        if (!skip_bw_average) {

            // Create channel: [ val(meta), [ bRep_bigwigs ] ]
            ch_binsize1
                .map { meta, bw ->
                    def meta_clone = meta.clone()
                    def antibody = meta.antibody ?: meta.input_control_of_antibody
                    meta_clone.id = meta_clone.id - ~/_bRep_.*$/
                    [ meta_clone.id, antibody, meta_clone, bw ]
                }
                .groupTuple(by: [0, 1])
                .map { id, antibody, metas, bws ->
                    def meta_clone = metas[0].clone()
                    meta_clone.averaged_brep = true
                    [meta_clone, bws.flatten()]
                }
                .set { ch_binsize1_brep_bw }

            //
            // MODULE: Average bigwigs across replicates
            //
            DEEPTOOLS_BIGWIGAVERAGE_BINSIZE1 (
                ch_binsize1_brep_bw,
                channel.value([[:], []])
            )
            ch_bw_avg_binsize1 = DEEPTOOLS_BIGWIGAVERAGE_BINSIZE1.out.bigwig
            ch_versions = ch_versions.mix(DEEPTOOLS_BIGWIGAVERAGE_BINSIZE1.out.versions.first())
        }

    }

    emit:
    bigwig_endo_rpm     = ch_bigwig_endo_rpm                           // channel: [ val(meta), [ bigwig ] ]
    bigwig_endo         = UCSC_BEDGRAPHTOBIGWIG_ENDO.out.bigwig    // channel: [ val(meta), [ bigwig ] ]
    bigwig_endo_rpm_lfc = ch_bigwig_endo_rpm_lfc                       // channel: [ val(meta), [ bigwig ] ]
    bigwig_avg_endo     = ch_bw_avg_endo        // channel: [ val(meta), [ bigwig ] ]
    bigwig_exo          = ch_bw_exo                                    // channel: [ val(meta), [ bigwig ] ]
    bigwig_binsize1     = ch_binsize1                                  // channel: [ val(meta), [ bigwig ] ]
    bigwig_binsize1_avg = ch_bw_avg_binsize1   // channel: [ val(meta), [ bigwig ] ]
    bedgraph_endo       = ch_bdg_all.filter { it -> it[0].genome == genome }      // channel: [ val(meta), [ bedgraph ] ]

    versions      = ch_versions                                     // channel: [ versions.yml ]
}