include { BAM_SPLIT_BY_STRAND                                       } from '../../../modules/local/bam_split_by_strand/main'
include { SAMTOOLS_INDEX                                            } from '../../../modules/nf-core/samtools/index/main'
include { BEDTOOLS_GENOMECOV                                        } from '../../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_MAKEWINDOWS                                      } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BIGTOOLS_BIGWIGAVERAGEOVERBED                             } from '../../../modules/local/bigtools/bigwigaverageoverbed/main'
include { BEDGRAPH_NORMALIZE                                            } from '../../../modules/local/bedgraph_normalize/main'
include { BEDGRAPH_SIGNAL_MINUS_INPUT                                   } from '../../../modules/local/bedgraph_signal_minus_input/main'
include { PARTITION_SMOOTH                                          } from '../../../modules/local/partition_smooth/main'
include { COLLECT_PARTITIONS                                          } from '../../../modules/local/collect_partitions/main'
include { BIGTOOLS_BEDGRAPHTOBIGWIG as BIGTOOLS_BEDGRAPHTOBIGWIG_WINDOWS } from '../../../modules/local/bigtools/bedgraphtobigwig/main'
include { BIGTOOLS_BEDGRAPHTOBIGWIG as BIGTOOLS_BEDGRAPHTOBIGWIG_PARTITIONS } from '../../../modules/local/bigtools/bedgraphtobigwig/main'
include { PARTITION_PLOT                                      } from '../../../modules/local/partition_plot/main'


workflow BAM_CREATE_SCAR_PARTITIONS {

    take:
    ch_bam                  // channel: [ val(meta), [ bam ] ]
    ch_chrom_sizes          // channel: [ bed ]
    ch_blacklist            // channel: [ val(meta), [ bed ] ]
    ch_initiation_zones     // channel: [ val(meta), [ bed ] ]
    rpm_use_flT2_total      // string: comma-separated list of antibodies for which to use flT2_total_mapped_reads instead of flT3_total_mapped_reads for RPM normalization


    main:

    ch_versions = Channel.empty()

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

    ch_bam = ch_f_bam.mix(ch_r_bam)

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

    //
    // MODULE: Calculate genome coverage
    //
    BEDTOOLS_GENOMECOV (
        ch_bam_scale,
        ch_chrom_sizes.map { it[1] },
        'bdg',
        true
    )
    ch_versions  = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    //
    // MODULE: Convert bedgraph to bigwig
    //
    BIGTOOLS_BEDGRAPHTOBIGWIG_WINDOWS (
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_chrom_sizes
    )
    ch_bigwig = BIGTOOLS_BEDGRAPHTOBIGWIG_WINDOWS.out.bigwig
    ch_versions = ch_versions.mix(BIGTOOLS_BEDGRAPHTOBIGWIG_WINDOWS.out.versions.first())

    //
    // MODULE: Create genomic windows
    //
    ch_windows = Channel.empty()
    BEDTOOLS_MAKEWINDOWS (
        ch_chrom_sizes
    )
    ch_windows = BEDTOOLS_MAKEWINDOWS.out.bed
    // count number of lines in the windows file
    ch_num_windows = ch_windows.map { meta, windows -> windows.countLines() }
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())

    //
    // MODULE: Calculate average coverage over windows
    //
    BIGTOOLS_BIGWIGAVERAGEOVERBED (
        ch_bigwig,
        ch_windows.first()
    )
    ch_bwaob = BIGTOOLS_BIGWIGAVERAGEOVERBED.out.bed
    ch_versions = ch_versions.mix(BIGTOOLS_BIGWIGAVERAGEOVERBED.out.versions.first())

    // TODO: print for debugging
    ch_bwaob
        .map { meta, bwaob ->
            "${meta}\t${bwaob}"
        }
        .collectFile( name: 'ch_bwaob.txt', newLine: true, sort: false, storeDir: "${params.outdir}/debug")


    // RPM normalization factors
    // num_windows is used to add a pseudocount to the RPM normalization factor (prevent division by zero in partition_smooth)
    ch_bwaob
        .combine(ch_num_windows)
        .map { meta, bwaob, num_windows ->
            def meta_clone = meta.clone()
            if (rpm_use_flT2_total && meta.antibody in rpm_use_flT2_total.split(',').collect { it.trim() }) {
                if (meta_clone.flT2_total_mapped_reads) {
                    meta_clone.norm_factor_val = 1e6 / (meta_clone.flT2_total_mapped_reads + num_windows)
                    meta_clone.norm_factor_val_used = 'flT2_total_mapped_reads'
                } else {
                    // Samples without spike-in wouldn't have flT2_total_mapped_reads, so we use flT1_total_mapped_reads instead
                    meta_clone.norm_factor_val = 1e6 / (meta_clone.flT1_total_mapped_reads + num_windows)
                    meta_clone.norm_factor_val_used = 'flT1_total_mapped_reads'
                }
            } else {
                meta_clone.norm_factor_val = 1e6 / (meta_clone.flT3_total_mapped_reads + num_windows)
                meta_clone.norm_factor_val_used = 'flT3_total_mapped_reads'
            }
            meta_clone.norm_factor_type = 'rpm'
            [ meta_clone, bwaob ]
        }
        .set { ch_bwaob_rpm }

    //
    // MODULE: Normalize strands
    //
    BEDGRAPH_NORMALIZE (
        ch_bwaob_rpm
    )
    ch_norm = BEDGRAPH_NORMALIZE.out.bedgraph
    ch_versions = ch_versions.mix(BEDGRAPH_NORMALIZE.out.versions.first())

    // for each of the strands, subtract the input from the sample
    ch_norm
        .branch { meta, bdg ->
            scar_with_input: meta.control
                return [ meta.control, meta.strand, meta, bdg ]
            input: !meta.control && meta.is_control
                return [ meta.id, meta.strand, bdg ]
        }
        .set { ch_norm_by_type }

    // create channel: [ val(meta), [ scar_bdg ], [ input_bdg ] ]
    ch_norm_by_type
        .scar_with_input
        .combine(ch_norm_by_type.input, by: [0, 1]) // combine by id and strand
        .map { input_id, strand, scar_meta, scar_bdg, input_bdg ->
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
        .collectFile( name: 'ch_norm_scar_input.txt', newLine: true, sort: false, storeDir: "${params.outdir}/debug") 


    //
    // MODULE: Substract input from sample
    //
    BEDGRAPH_SIGNAL_MINUS_INPUT (
        ch_norm_scar_input
    )
    ch_bdg_smi = BEDGRAPH_SIGNAL_MINUS_INPUT.out.bedgraph
    ch_versions = ch_versions.mix(BEDGRAPH_SIGNAL_MINUS_INPUT.out.versions.first())

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
        .collectFile( name: 'ch_norm_and_smi.txt', newLine: true, sort: false, storeDir: "${params.outdir}/debug")

    //
    // MODULE: Calculate partitions (RFD)
    //
    PARTITION_SMOOTH (
        ch_norm_and_smi,
        params.scar_radius,
        params.scar_dradius,
        params.scar_zradius
    )
    ch_rfd = PARTITION_SMOOTH.out.rfd
    ch_versions = ch_versions.mix(PARTITION_SMOOTH.out.versions.first())

    // TODO: print for debugging
    ch_rfd
        .map { meta, rfd ->
            "${meta}\t${rfd}"
        }
        .collectFile( name: 'ch_rfd.txt', newLine: true, sort: false, storeDir: "${params.outdir}/debug")


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
            meta_clone.removeAll { it.key in ['norm_factor_val', 'norm_factor_val_used', 'norm_factor_type', 'signal_minus_input'] }
            [ meta_clone, meta, norm_or_smi_fwd, norm_or_smi_rev ]
        }
        .set { ch_norm_and_smi_to_combine }

    // TODO: print for debugging
    ch_norm_and_smi_to_combine
        .map { meta, meta_norm_or_smi, norm_or_smi_fwd, norm_or_smi_rev ->
            "${meta}\t${meta_norm_or_smi}\t${norm_or_smi_fwd}\t${norm_or_smi_rev}"
        }
        .collectFile( name: 'ch_norm_and_smi_to_combine.txt', newLine: true, sort: false, storeDir: "${params.outdir}/debug")

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
        .set { ch_partitions }

    // TODO: print for debugging
    ch_partitions
        .map { meta, windows, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ->
            "${meta}\t${windows}\t${bwaob_fwd}\t${bwaob_rev}\t${norm_or_smi_fwd}\t${norm_or_smi_rev}\t${rfd}"
        }
        .collectFile( name: 'ch_partitions.txt', newLine: true, sort: false, storeDir: "${params.outdir}/debug")

    //
    // MODULE: Collect partitions
    //
    COLLECT_PARTITIONS (
        ch_partitions
    )
    ch_versions = ch_versions.mix(COLLECT_PARTITIONS.out.versions.first())

    COLLECT_PARTITIONS
        .out
        .tsv
        .branch { meta, tsv ->
            scar_with_input: !meta.is_control && !meta.signal_minus_input
                return [ meta.control, meta, tsv ]
            input: meta.is_control
                return [ meta.id, tsv ]
            minusinput: meta.signal_minus_input
                return [ meta.id, tsv ]
        }
        .set { ch_partitions_by_type }

    ch_partitions_by_type
        .scar_with_input
        .combine(ch_partitions_by_type.input, by: 0) // this creates channel: [ input_id, meta_scar, scar_tsv, input_tsv ]
        .map { input_id, meta_scar, scar_tsv, input_tsv ->
            [ meta_scar.id, meta_scar, scar_tsv, input_tsv ]
        }
        .combine(ch_partitions_by_type.minusinput, by: 0) // this creates channel: [ scar_id, meta_scar, scar_tsv, input_tsv, minusinput_tsv ]
        .map { scar_id, meta_scar, scar_tsv, input_tsv, minusinput_tsv ->
            def okseq = meta_scar.okseq_part_file ? file(meta_scar.okseq_part_file) : null
            [ scar_id, meta_scar, scar_tsv, input_tsv, minusinput_tsv, okseq ]
        }
        .set { ch_partitions_to_plot }

    // TODO: print for debugging
    ch_partitions_to_plot
        .map { scar_id, meta_scar, scar_tsv, input_tsv, minusinput_tsv, okseq ->
            "${scar_id}\t${meta_scar}\t${scar_tsv}\t${input_tsv}\t${minusinput_tsv}\t${okseq}"
        }
        .collectFile( name: 'ch_partitions_to_plot.txt', newLine: true, sort: false, storeDir: "${params.outdir}/debug")


    //
    // MODULE: Plot the final partition
    //
    PARTITION_PLOT (
        ch_partitions_to_plot,
        ch_blacklist,
        ch_initiation_zones
    )
    ch_versions = ch_versions.mix(PARTITION_PLOT.out.versions.first())

    //
    // MODULE: Convert the final partition bedgraph to bigwig
    //
    BIGTOOLS_BEDGRAPHTOBIGWIG_PARTITIONS (
        COLLECT_PARTITIONS.out.bdg,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(BIGTOOLS_BEDGRAPHTOBIGWIG_PARTITIONS.out.versions.first())

    emit:
    tab      = PARTITION_SMOOTH.out.rfd       // channel: [ val(meta), [ tab ] ]
    versions = ch_versions                    // channel: [ versions.yml ]
}

