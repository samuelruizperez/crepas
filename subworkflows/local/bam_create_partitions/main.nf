include { BAM_SPLIT_BY_STRAND                                       } from '../../../modules/local/bam_split_by_strand/main'
include { SAMTOOLS_INDEX                                            } from '../../../modules/nf-core/samtools/index/main'
include { BEDTOOLS_GENOMECOV                                        } from '../../../modules/nf-core/bedtools/genomecov/main'
include { FILE_SORT as BEDGRAPH_SORT                                } from '../../../modules/local/file_sort/main'
include { BEDTOOLS_MAKEWINDOWS                                      } from '../../../modules/nf-core/bedtools/makewindows/main'
include { FILE_SORT as WINDOWS_SORT                                 } from '../../../modules/local/file_sort/main'
include { UCSC_BIGWIGAVERAGEOVERBED                                 } from '../../../modules/nf-core/ucsc/bigwigaverageoverbed/main'
include { FILE_SORT as BWAOB_SORT                                   } from '../../../modules/local/file_sort/main'
include { BEDGRAPH_NORMALIZE as BWAOB_NORMALIZE                     } from '../../../modules/local/bedgraph_normalize/main'
include { BEDGRAPH_SIGNAL_MINUS_INPUT                               } from '../../../modules/local/bedgraph_signal_minus_input/main'
include { PARTITION_OR_RFD_SMOOTH                                   } from '../../../modules/local/partition_or_rfd_smooth/main'
include { COLLECT_PARTITIONS                                        } from '../../../modules/local/collect_partitions/main'
include { PARTITION_AVERAGE                                         } from '../../../modules/local/partition_average/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_WINDOWS    } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_PARTITIONS } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { RFD_TO_IZ                                               } from '../../../modules/local/rfd_to_iz/main'
include { PARTITION_PLOT                                            } from '../../../modules/local/partition_plot/main'


workflow BAM_CREATE_PARTITIONS {

    take:
    ch_bam                  // channel: [ val(meta), [ bam ] ]
    ch_chrom_sizes          // channel: [ bed ]
    ch_blacklist            // channel: [ val(meta), [ bed ] ]
    ch_okseq_rfd_file       // channel: [ val(meta), [ bed ] ]
    ch_initiation_zones     // channel: [ val(meta), [ bed ] ]
    smooth_radius
    derivative_radius
    zero_crossing_radius

    main:

    ch_versions = channel.empty()

    // TODO: print for debugging
    ch_bam
        .map { meta, bam ->
            "${meta}\t${bam}"
        }
        .collectFile( name: '1_scar_ch_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    //
    // MODULE: Split BAMs by strand (forward and reverse)
    //
    BAM_SPLIT_BY_STRAND ( ch_bam )
    ch_versions = ch_versions.mix(BAM_SPLIT_BY_STRAND.out.versions.first())

    // Add strand to the meta information
    BAM_SPLIT_BY_STRAND
        .out
        .f_bam
        .map {
            meta, f_bam ->
                def meta_clone = meta.clone()
                meta_clone.strand = 'forward'
                [ meta_clone, f_bam ]
        }
        .set { ch_f_bam }

    // TODO: print for debugging
    ch_f_bam
        .map { meta, f_bam ->
            "${meta}\t${f_bam}"
        }
        .collectFile( name: '2_scar_ch_f_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    BAM_SPLIT_BY_STRAND
        .out
        .r_bam
        .map {
            meta, r_bam ->
                def meta_clone = meta.clone()
                meta_clone.strand = 'reverse'
                [ meta_clone, r_bam ]
        }
        .set { ch_r_bam }
    
    // TODO: print for debugging
    ch_r_bam
        .map { meta, r_bam ->
            "${meta}\t${r_bam}"
        }
        .collectFile( name: '3_scar_ch_r_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    ch_bam = ch_f_bam.mix(ch_r_bam)

    // TODO: print for debugging
    ch_bam
        .map { meta, bam ->
            "${meta}\t${bam}"
        }
        .collectFile( name: '4_scar_ch_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Index BAM files per strand
    //
    SAMTOOLS_INDEX (
        ch_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Creating channel: [ val(meta), [ bam ], [ scale ] ] 
    ch_bam
        .map {
            meta, bam ->
                [ meta, bam, 1 ]
        }
        .set { ch_bam_scale }

    // TODO: print for debugging
    ch_bam_scale
        .map { meta, bam, scale ->
            "${meta}\t${bam}\t${scale}"
        }
        .collectFile( name: '5_scar_ch_bam_scale.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Calculate genome coverage
    //
    BEDTOOLS_GENOMECOV (
        ch_bam_scale,
        ch_chrom_sizes.map { it -> it[1] },
        'bdg',
        true
    )
    ch_versions  = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    //
    // MODULE: Sort the bedgraph so that it works with bedgraphtobigwig
    //
    BEDGRAPH_SORT (
        BEDTOOLS_GENOMECOV.out.genomecov,
        'bedGraph'
    )
    ch_versions = ch_versions.mix(BEDGRAPH_SORT.out.versions.first())

    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_WINDOWS (
        BEDGRAPH_SORT.out.sorted,
        ch_chrom_sizes.map { it -> it[1] }
    )
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG_WINDOWS.out.bigwig
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_WINDOWS.out.versions.first())

    // TODO: print for debugging
    ch_bigwig
        .map { meta, bigwig ->
            "${meta}\t${bigwig}"
        }
        .collectFile( name: '6_scar_ch_bigwig.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Create genomic windows
    //
    ch_windows = channel.empty()
    BEDTOOLS_MAKEWINDOWS (
        ch_chrom_sizes
    )
    ch_windows = BEDTOOLS_MAKEWINDOWS.out.bed

    //
    // MODULE: Sort windows
    //
    WINDOWS_SORT (
        ch_windows,
        'bed'
    )
    ch_windows = WINDOWS_SORT.out.sorted

    // count number of lines in the windows file
    ch_num_windows = ch_windows.map { meta, windows -> windows.countLines() }
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

    // TODO: print for debugging
    ch_num_windows
        .map { num_windows ->
            "${num_windows}"
        }
        .collectFile( name: '7_scar_ch_num_windows.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    ch_bigwig
        .combine(ch_windows)
        .map { meta, bigwig, meta_windows, windows ->
            [ meta, windows, bigwig ]
        }
        .set { ch_windows_bigwig }

    //
    // MODULE: Calculate average coverage over windows
    //
    UCSC_BIGWIGAVERAGEOVERBED (
        ch_windows_bigwig
    )
    ch_bwaob = UCSC_BIGWIGAVERAGEOVERBED.out.tab
    ch_versions = ch_versions.mix(UCSC_BIGWIGAVERAGEOVERBED.out.versions.first())

    // TODO: print for debugging
    ch_bwaob
        .map { meta, bwaob ->
            "${meta}\t${bwaob}"
        }
        .collectFile( name: '8_scar_ch_bwaob.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Sort BWAOB
    //
    BWAOB_SORT (
        ch_bwaob,
        'tab'
    )
    ch_bwaob = BWAOB_SORT.out.sorted
    ch_versions = ch_versions.mix(BWAOB_SORT.out.versions.first())



    // RPM normalization factors
    // num_windows is used to add a pseudocount to the RPM normalization factor (prevent division by zero in partition_or_rfd_smooth)
    ch_bwaob
        .combine(ch_num_windows)
        .map { meta, bwaob, num_windows ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1e6 / (meta[meta.ref_total_mapped_reads_for_rpm] + num_windows)
            meta_clone.norm_factor_type = 'rpm'
            [ meta_clone, bwaob ]
        }
        .set { ch_bwaob_rpm }

    // TODO: print for debugging
    ch_bwaob_rpm
        .map { meta, bwaob ->
            "${meta}\t${bwaob}"
        }
        .collectFile( name: '9_scar_ch_bwaob_rpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Normalize the 4th column of the bwaob file (sum of values over all bases covered)
    //
    BWAOB_NORMALIZE (
        ch_bwaob_rpm,
        'tab'
    )
    ch_norm = BWAOB_NORMALIZE.out.normalized
    ch_versions = ch_versions.mix(BWAOB_NORMALIZE.out.versions.first())

    // TODO: print for debugging
    ch_norm
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: '10_scar_ch_norm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    // for each of the strands, subtract the input from the sample
    ch_norm
        .branch { meta, bdg ->
            scar_with_ipcontrol: meta.input_control
                return [ meta.input_control, meta.antibody, meta.strand, meta, bdg ]
            ipcontrol: !meta.input_control && meta.is_input_control
                return [ meta.id, meta.input_control_of_antibody, meta.strand, bdg ]
        }
        .set { ch_norm_by_type }

    // create channel: [ val(meta), [ scar_bdg ], [ input_bdg ] ]
    ch_norm_by_type
        .scar_with_ipcontrol
        .combine(ch_norm_by_type.ipcontrol, by: [0, 1, 2]) // combine by id, antibody, and strand
        .map { ipcontrol_id, antibody, strand, scar_meta, scar_bdg, input_bdg ->
            def meta_clone = scar_meta.clone()
                meta_clone.signal_minus_input = true
                [ meta_clone, scar_bdg, input_bdg ]
        }
        .set { ch_norm_scar_input }


    // TODO: print for debugging
    ch_norm_scar_input
        .map { meta, scar_bdg, input_bdg ->
            "${meta}\t${scar_bdg}\t${input_bdg}"
        }
        .collectFile( name: '11_scar_ch_norm_scar_input.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    //
    // MODULE: Substract input from sample
    //
    BEDGRAPH_SIGNAL_MINUS_INPUT (
        ch_norm_scar_input
    )
    ch_bdg_smi = BEDGRAPH_SIGNAL_MINUS_INPUT.out.bedgraph
    ch_versions = ch_versions.mix(BEDGRAPH_SIGNAL_MINUS_INPUT.out.versions.first())

    // TODO: print for debugging
    ch_bdg_smi
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: '12_scar_ch_bdg_smi.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    // create channel: [ val(meta), [ bdg_fwd ], [ bdg_rev ] ]
    ch_norm
        .mix(ch_bdg_smi)
        // copy meta and remove meta.strand to then merge fwd and rev by meta
        .map {
            meta, bdg ->
                def meta_clone = meta.clone()
                meta_clone.remove('strand')
                [ meta_clone, meta, bdg ]
        }
        .branch { meta_clone, meta, bdg ->
            forward: meta.strand == 'forward'
            reverse: meta.strand == 'reverse'
        }
        .set { ch_norm_and_smi }

    ch_norm_and_smi.forward
        .combine(ch_norm_and_smi.reverse, by: 0)
        .map { meta_clone, meta_fwd, bdg_fwd, meta_rev, bdg_rev ->
            [ meta_clone, bdg_fwd, bdg_rev ]
        }
        .set { ch_norm_and_smi }

    // TODO: print for debugging
    ch_norm_and_smi
        .map { meta, bdg_fwd, bdg_rev ->
            "${meta}\t${bdg_fwd}\t${bdg_rev}"
        }
        .collectFile( name: '13_scar_ch_norm_and_smi.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    // Create channel: [ val(meta), val(partition_or_rfd), [ f_tab ], [ r_tab ] ]
    ch_norm_and_smi
        .map { meta, bdg_fwd, bdg_rev ->
            // 'partition' is for SCAR-seq and 'RFD' for OK-seq
            def partition_or_rfd = meta.exp_type == 'SCAR-seq' ? 'partition' : meta.exp_type == 'OK-seq' ? 'RFD' : null
            [ meta, partition_or_rfd, bdg_fwd, bdg_rev ]
        }
        .set { ch_part_norm_and_smi }

    //
    // MODULE: Calculate partitions (RFD)
    //
    PARTITION_OR_RFD_SMOOTH (
        ch_part_norm_and_smi,
        smooth_radius,
        derivative_radius,
        zero_crossing_radius
    )
    ch_rfd = PARTITION_OR_RFD_SMOOTH.out.rfd
    ch_versions = ch_versions.mix(PARTITION_OR_RFD_SMOOTH.out.versions.first())

    // TODO: print for debugging
    ch_rfd
        .map { meta, rfd ->
            "${meta}\t${rfd}"
        }
        .collectFile( name: '14_scar_ch_rfd.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    // Prepare bwaob channel for combine()
    ch_bwaob
        .map { meta, bwaob ->
            def meta_clone = meta.clone()
            meta_clone.remove('strand')
            [ meta.strand, meta_clone, bwaob ]
        }
        .branch { strand, meta, bwaob ->
            forward: strand == 'forward'
                return [ meta, bwaob ]
            reverse: strand == 'reverse'
                return [ meta, bwaob ]
        }
        .set { ch_bwaob_strands }

    // Prepare norm and smi channel for combine()
    ch_norm_and_smi
        .map { meta, norm_or_smi_fwd, norm_or_smi_rev ->
            def meta_clone = meta.clone()
            meta_clone.removeAll { it -> it.key in ['norm_factor_val', 'norm_factor_type', 'signal_minus_input'] }
            [ meta_clone, meta, norm_or_smi_fwd, norm_or_smi_rev ]
        }
        .set { ch_norm_and_smi_to_combine }


    // TODO: print for debugging
    ch_norm_and_smi_to_combine
        .map { meta, meta_norm_or_smi, norm_or_smi_fwd, norm_or_smi_rev ->
            "${meta}\t${meta_norm_or_smi}\t${norm_or_smi_fwd}\t${norm_or_smi_rev}"
        }
        .collectFile( name: '15_scar_ch_norm_and_smi_to_combine.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    // Create channel: [ meta, windows, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ]
    ch_bwaob_strands.forward
        .combine(ch_bwaob_strands.reverse, by: 0) // this creates channel: [ meta, bwaob_fwd, bwaob_rev ]
        .combine(ch_norm_and_smi_to_combine, by: 0) // this creates channel: [ meta, bwaob_fwd, bwaob_rev, meta_norm_or_smi, norm_or_smi_fwd, norm_or_smi_rev ]
        .map { meta, bwaob_fwd, bwaob_rev, meta_norm_or_smi, norm_or_smi_fwd, norm_or_smi_rev ->
            [ meta_norm_or_smi, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev ]
        }
        .combine(ch_rfd, by: 0) // this creates channel: [ meta, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ]
        .combine(ch_windows) // this creates channel: [ meta, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd, windows ]
        .map { meta, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd, meta_windows, windows ->
            [ meta, windows, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ]
        }
        .set { ch_to_collect }

    // TODO: print for debugging
    ch_to_collect
        .map { meta, windows, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ->
            "${meta}\t${windows}\t${bwaob_fwd}\t${bwaob_rev}\t${norm_or_smi_fwd}\t${norm_or_smi_rev}\t${rfd}"
        }
        .collectFile( name: '16_scar_ch_partitions.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Collect partitions
    //
    COLLECT_PARTITIONS (
        ch_to_collect
    )
    ch_partitions = COLLECT_PARTITIONS.out.tsv
    ch_versions = ch_versions.mix(COLLECT_PARTITIONS.out.versions.first())


    // Create channel: [ val(meta), partitions_brep ]
    ch_partitions_brep = channel.empty()
    ch_partitions
        .map { meta, partition ->
            def meta_clone = meta.clone()
            def antibody = meta.antibody ?: meta.input_control_of_antibody
            meta_clone.id = meta.id - ~/_bRep_.*$/
            meta_clone.input_control = meta.input_control - ~/_bRep_.*$/
            [ meta_clone.id, antibody, meta.signal_minus_input, meta_clone, partition ]
        }
        .groupTuple(by: [0, 1, 2])
        // remove elements where there is only one biological replicate
        .filter { id, antibody, smi, metas, partitions ->
            partitions.size() > 1
        }
        .map { id, antibody, smi, metas, partitions ->
            def meta_clone = metas[0].clone()
            meta_clone.averaged_brep = true
            [meta_clone, partitions.flatten()]
        }
        .set { ch_partitions_brep }

    // TODO: print for debugging
    ch_partitions_brep
        .map { meta, partitions ->
            "${meta}\t${partitions}"
        }
        .collectFile( name: '18_scar_ch_partitions_brep.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Create average partition across biological replicates
    //
    PARTITION_AVERAGE (
        ch_partitions_brep
    )
    ch_versions = ch_versions.mix(PARTITION_AVERAGE.out.versions.first())


    ch_partitions
        .mix(PARTITION_AVERAGE.out.tsv)
        .filter { it -> it[0].exp_type == 'OK-seq' }
        .set { ch_okseq }
    
    //
    // MODULE: Process OK-seq RFD file to get initiation zones
    //
    // Note: this IZ file is not used for plotting in the PARTITION_PLOT module,
    // it is only created here when an OK-seq sample is being processed
    // The one used for plotting is generated in PREPARE_GENOME
    RFD_TO_IZ (
        ch_okseq,
        ch_blacklist,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(RFD_TO_IZ.out.versions)

    ch_partitions
        .mix(PARTITION_AVERAGE.out.tsv)
        .branch { meta, tsv ->
            scar_with_ipcontrol: !meta.is_input_control && !meta.signal_minus_input
                return [ meta.input_control, meta, tsv ]
            ipcontrol: meta.is_input_control
                return [ meta.id, tsv ]
            minusipcontrol: meta.signal_minus_input
                return [ meta.id, tsv ]
        }
        .set { ch_partitions_by_type }

    ch_partitions_by_type
        .scar_with_ipcontrol
        .combine(ch_partitions_by_type.ipcontrol, by: 0)
        .map { ipcontrol_id, meta_scar, scar_tsv, input_tsv ->
            [ meta_scar.id, meta_scar, scar_tsv, input_tsv ]
        }
        .combine(ch_partitions_by_type.minusipcontrol, by: 0)
        .map { scar_id, meta_scar, scar_tsv, input_tsv, minusinput_tsv ->
            [ meta_scar, scar_tsv, input_tsv, minusinput_tsv ]
        }
        .set { ch_partitions_to_plot }


    // TODO: print for debugging
    ch_partitions_to_plot
        .map { meta_scar, scar_tsv, input_tsv, minusinput_tsv ->
            "${meta_scar}\t${scar_tsv}\t${input_tsv}\t${minusinput_tsv}"
        }
        .collectFile( name: '17_scar_ch_partitions_to_plot.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Plot the final partition
    //
    PARTITION_PLOT (
        ch_partitions_to_plot,
        ch_blacklist,
        ch_okseq_rfd_file,
        ch_initiation_zones,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(PARTITION_PLOT.out.versions.first())

    //
    // MODULE: Convert the final partition bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_PARTITIONS (
        COLLECT_PARTITIONS.out.bdg.mix(PARTITION_AVERAGE.out.bdg),
        ch_chrom_sizes.map { it -> it[1] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_PARTITIONS.out.versions.first())

    emit:
    tab      = PARTITION_OR_RFD_SMOOTH.out.rfd       // channel: [ val(meta), [ tab ] ]
    versions = ch_versions                    // channel: [ versions.yml ]
}

