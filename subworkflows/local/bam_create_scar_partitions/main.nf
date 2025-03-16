include { BAM_SPLIT_BY_STRAND                                       } from '../../../modules/local/bam_split_by_strand/main'
include { SAMTOOLS_INDEX                                            } from '../../../modules/nf-core/samtools/index/main'
include { BEDTOOLS_GENOMECOV                                        } from '../../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_SLOP                                             } from '../../../modules/nf-core/bedtools/slop/main'
include { UCSC_BEDCLIP                                              } from '../../../modules/nf-core/ucsc/bedclip/main'
include { BEDTOOLS_MAKEWINDOWS                                      } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BED_SPLIT_BY_CHROMOSOME                                   } from '../../../modules/local/bed_split_by_chromosome/main'
include { UCSC_BIGWIGAVERAGEOVERBED                                 } from '../../../modules/nf-core/ucsc/bigwigaverageoverbed/main'
include { CPM_CALCULATION as CPM_CALCULATION_SAMPLES                } from '../../../modules/local/cpm_calculation/main'
include { CPM_CALCULATION as CPM_CALCULATION_INPUTS                 } from '../../../modules/local/cpm_calculation/main'
include { NORMALIZE_STRANDS                                         } from '../../../modules/local/normalize_strands/main'
include { SUBSTRACT_INPUT                                           } from '../../../modules/local/substract_input/main'
include { PARTITION_SMOOTH                                          } from '../../../modules/local/partition_smooth/main'
include { COLLECT_PARTITIONS_BY_CHROMOSOME                          } from '../../../modules/local/collect_partitions_by_chromosome/main'
include { FINAL_PARTITION_BEDGRAPH                                  } from '../../../modules/local/final_partition_bedgraph/main'
include { FILE_SORT as FILE_SORT_WINDOWS                            } from '../../../modules/local/file_sort/main'
include { FILE_SORT as FILE_SORT_PARTITIONS                         } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_WINDOWS    } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_PARTITIONS } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { FINAL_PARTITION_PLOT                                      } from '../../../modules/local/final_partition_plot/main'


workflow BAM_CREATE_SCAR_PARTITIONS {

    take:
    ch_bam                  // channel: [ val(meta), [ bam ] ]
    ch_chrom_sizes          // channel: [ bed ]
    ch_blacklist            // channel: [ val(meta), [ bed ] ]
    ch_initiation_zones     // channel: [ val(meta), [ bed ] ]
    //ch_scaffolds          // channel: [scaffolds]

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
                [ meta + ['strand':'forward'], f_bam ]
        }
        .set { ch_f_bam }

    BAM_SPLIT_BY_STRAND
        .out
        .r_bam
        .map {
            meta, r_bam ->
                [ meta + ['strand':'reverse'], r_bam ]
        }
        .set { ch_r_bam }

    ch_bam = ch_f_bam.mix(ch_r_bam)

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
    // MODULE: Increase the size of each feature
    //
    BEDTOOLS_SLOP (
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_chrom_sizes.map { it[1] }
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SLOP.out.versions.first())

    //
    // MODULE: Remove records that are out of annotated chromosome ranges
    //
    UCSC_BEDCLIP (
        BEDTOOLS_SLOP.out.bed,
        ch_chrom_sizes.map { it[1] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions.first())

    //
    // MODULE: Sort the bedgraph file
    //
    FILE_SORT_WINDOWS (
        UCSC_BEDCLIP.out.bedgraph,
        'sorted'
    )
    ch_versions = ch_versions.mix(FILE_SORT_WINDOWS.out.versions.first())

    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_WINDOWS (
        FILE_SORT_WINDOWS.out.sorted,
        ch_chrom_sizes.map { it[1] }
    )
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG_WINDOWS.out.bigwig
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_WINDOWS.out.versions.first())

    //
    // MODULE: Create chromosome windows
    //

    // TODO: windows are created and split even when not needed (no scarseq samples)
    // this is an ugly workaround
    // https://github.com/nextflow-io/nextflow/discussions/5102#discussioncomment-9939140
    // def make_windows = false
    // //ch_chrom_sizes_endo_ss = ch_chrom_sizes_endo
    // ch_bigwig.count().map { n ->
    //     if (n > 0) {
    //         make_windows = true
    //     }
    // }

    ch_windows = Channel.empty()
    //if (make_windows) {
        BEDTOOLS_MAKEWINDOWS (
            ch_chrom_sizes
        )
        ch_windows = BEDTOOLS_MAKEWINDOWS.out.bed
        ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())
    //}

    // Create a channel with each chromosome to iterate over
    ch_chroms = Channel.empty()
    ch_chrom_sizes
        .map {
            meta, bed ->
                bed.splitCsv(header:false, sep:'\t')
        }
        .flatMap { chrom_list ->
        chrom_list.collect { it[0] }
        }
        .set { ch_chroms }


    //
    // MODULE: Split BED (windows) by chromosome
    //
    BED_SPLIT_BY_CHROMOSOME (
        ch_chroms,
        ch_windows.first()
    )
    ch_chroms   = BED_SPLIT_BY_CHROMOSOME.out.bed
    ch_versions = ch_versions.mix(BED_SPLIT_BY_CHROMOSOME.out.versions.first())


    // Create channel by combining the bigwig files with the chromosomes (all chroms per bigwig)
    ch_bigwig
        .combine(ch_chroms)
        // Append the chromosome meta information to the bigwig meta information
        .map {
            meta, bigWig, chr_meta, chr_bed ->
                [ meta + chr_meta, bigWig, chr_bed ]
        }
        .set { ch_bigwig_chroms }
    
    // Separate the channel again because UCSC_BIGWIGAVERAGEOVERBED
    // expects a channel of bigwig files and a channel of chromosome beds
    ch_bigwig_chroms
        .map {
            meta, bigwig, chrom ->
                [ meta, bigwig ]
        }
        .set {ch_bw_combs}

    // Add bigwig meta information to the chromosome bed files
    ch_bigwig_chroms
        .map {
            meta, bigwig, chrom ->
                [ meta, chrom ]
            }
        .set { ch_chroms_combs }


    //
    // MODULE: Calculate average coverage over windows
    //
    UCSC_BIGWIGAVERAGEOVERBED (
        ch_chroms_combs,
        ch_bw_combs.map{ it[1] }
    )
    ch_bwaob = UCSC_BIGWIGAVERAGEOVERBED.out.tab
    ch_versions = ch_versions.mix(UCSC_BIGWIGAVERAGEOVERBED.out.versions.first())


    // Separate samples and input controls
    ch_bwaob
        .map {
            meta, chroms ->
                meta.control ? null : [ meta, chroms ]
        }
        .set { ch_inputs }

    ch_bwaob
        .map {
            meta, chroms ->
                meta.control ? [ meta, chroms ] : null
        }
        .set { ch_samples }

    // for each ch_samples.meta.id, concatenate the files (use collectFile)
    // remove strand information from the meta information
    ch_inputs
        .map {
            meta, chroms ->
                def meta_clone = meta.clone()
                meta_clone.remove('strand')
                [ meta_clone ]
        }
        .unique()
        .map { it[0] }
        .map { it -> [ it.id, it ] }
        .set { ch_i_uniq_meta }


    ch_inputs
        .collectFile(newLine: false, sort: true) { //, storeDir: "${params.outdir}/${params.aligner}/mergedLibrary/scarseq/windows_cat") {
            meta, tabs ->
                // to do it per id and strand:
                // [ "${meta.id}.${meta.strand}.inputs_windows.tab", tabs ]
                // to do it per id only
                [ "${meta.id}.bwaob.windows.tab", tabs ]
        }
        // readd meta information based on meta.id in the filename
        .map {
            tab ->
                def id = tab.name.split("\\.")[0]
                [ id, tab ]
        }
        .join(ch_i_uniq_meta, by: 0)
        // flip meta and tab to have meta first
        .map {
            id, tab, meta ->
                [ meta, tab ]
        }
        .set { ch_inputs_bed }


    ch_samples
        .map {
            meta, chroms ->
                def meta_clone = meta.clone()
                meta_clone.remove('strand')
                [ meta_clone ]
        }
        .unique()
        .map { it[0] }
        .map { it -> [ it.id, it ] }
        .set { ch_s_uniq_meta }



    ch_samples
        .collectFile(newLine: false, sort: true) { //, storeDir: "${params.outdir}/${params.aligner}/mergedLibrary/scarseq/windows_cat") {
            meta, tabs ->
                // to do it per id and strand:
                // [ "${meta.id}.${meta.strand}.samples_windows.tab", tabs ]
                // to do it per id only
                [ "${meta.id}.bwaob.windows.tab", tabs ]
        }
        // re-add meta information based on meta.id in the filename
        .map {
            tab ->
                def id = tab.name.split("\\.")[0]
                [ id, tab ]
        }
        .join(ch_s_uniq_meta, by: 0)
        // flip meta and tab to have meta first
        .map {
            id, tab, meta ->
                [ meta, tab ]
        }
        .set { ch_samples_bed }
    

    //
    // MODULE: Calculate CPM
    //
    CPM_CALCULATION_SAMPLES (
        ch_samples_bed
    )
    ch_versions = ch_versions.mix(CPM_CALCULATION_SAMPLES.out.versions.first())

    ch_cpm_samples = CPM_CALCULATION_SAMPLES.out.cpm
        .map {
            meta, cpm ->
                [ meta.id, cpm.splitCsv(header:false)[0][0] ]
        }

    CPM_CALCULATION_INPUTS (
        ch_inputs_bed
    )
    ch_versions = ch_versions.mix(CPM_CALCULATION_INPUTS.out.versions.first())

    ch_cpm_inputs = CPM_CALCULATION_INPUTS.out.cpm
        .map {
            meta, cpm ->
                [ meta.id, cpm.splitCsv(header:false)[0][0] ] // TODO: check if flatten is correct
        }

    // concat samples and inputs and add their corresponding cpm
    ch_samples
        .map {
            meta, chroms ->
                [ meta.id, meta, chroms ]
        }
        .combine(ch_cpm_samples, by: 0)
        .concat(
            ch_inputs
            .map {
                meta, chroms ->
                    [ meta.id, meta, chroms ]
                }
            .combine(ch_cpm_inputs, by: 0)
        )
        .map {
            id, meta, chroms, cpm ->
                [ meta, chroms, cpm ]
        }
        .set { ch_samples_inputs_cpm }

    //
    // MODULE: Normalize strands
    //
    NORMALIZE_STRANDS (
        ch_samples_inputs_cpm,
    )
    ch_normalized_strands = NORMALIZE_STRANDS.out.tab
    ch_versions = ch_versions.mix(NORMALIZE_STRANDS.out.versions.first())

    // for each of the chromosomes, and each of the strands, subtract the input from the sample
    //ch_norm_s_i_cpm = Channel.empty()

    // Step 1: Separate samples and controls
    ch_normalized_strands
        .filter { meta, tab, cpm -> meta.control } // this is the groovy way to check if meta.control is empty/not null
        .set { ch_samples }

    ch_normalized_strands
        .filter { meta, tab, cpm -> !meta.control }
        .set { ch_controls }


    // Step 2: Combine samples and controls based on matching chromosome and strand
    ch_samples
        .combine(ch_controls)
        .filter { sampleMeta, sampleTab, sampleCpm, controlMeta, controlTab, controlCpm ->
            sampleMeta.chr == controlMeta.chr &&
            sampleMeta.strand == controlMeta.strand &&
            sampleMeta.control == controlMeta.id
        }
        .map { sampleMeta, sampleTab, sampleCpm, controlMeta, controlTab, controlCpm ->
            tuple(sampleMeta, sampleTab, controlTab, sampleCpm)
        }
        .set { ch_norm_s_i_cpm }


    //
    // MODULE: Substract input from sample
    //
    SUBSTRACT_INPUT (
        ch_norm_s_i_cpm
    )
    ch_substracted = SUBSTRACT_INPUT.out.tab
    ch_versions = ch_versions.mix(SUBSTRACT_INPUT.out.versions.first())


    // concat normalized output and substracted output
    ch_normalized_strands
        .map {
            meta, tab, cpm ->
                [ meta, tab ]
        }
        .concat(ch_substracted)
        // copy meta and remove meta.strand from meta clone
        .map {
            meta, tab ->
                def meta_clone = meta.clone()
                meta_clone.remove('strand')
                [ meta_clone, meta, tab ]
        }
        .branch { meta_clone, meta, tab ->
            forward: meta.strand == 'forward'
            reverse: meta.strand == 'reverse'
        }
        .set { ch_norm_and_subs }  // Assign to new channel

    ch_norm_and_subs.forward
        .combine(ch_norm_and_subs.reverse, by: 0)
        .map { meta_clone, meta1, tab1, meta2, tab2 ->
            [ meta_clone, tab1, tab2 ]
        }
        .set { ch_norm_and_subs }


    //
    // MODULE: Smooth the partition
    //
    PARTITION_SMOOTH (
        ch_norm_and_subs,
        params.scar_radius,
        params.scar_dradius,
        params.scar_zradius
    )
    ch_part_smooth = PARTITION_SMOOTH.out.rfd
    ch_versions = ch_versions.mix(PARTITION_SMOOTH.out.versions.first())


    ch_bwaob
        .concat(ch_normalized_strands.map { it -> [ it[0], it[1] ] }) // remove cpm
        .concat(ch_part_smooth.filter { meta, tab -> !meta.minusinput }) // remove minusinput
        .map { meta, tab -> [ meta.id, meta.chr, [meta, tab] ] }
        .groupTuple( by: [0, 1] )
        .map { id, chr, group ->
            def sortedGroup = group
                .sort { a, b ->
                    // Sorting criteria:
                    a[0].RFD <=> b[0].RFD ?:           // Sort by RFD (true vs false)
                    a[0].cpm <=> b[0].cpm ?:          // Sort by cpm (true vs false)
                    a[0].strand <=> b[0].strand    // Sort by strand (forward vs reverse)
                }

            def primaryMeta = sortedGroup[0][0]      // Take the first meta (after sorting)
            def filesList = sortedGroup*.get(1)      // Extract list of files

            [id, chr, [primaryMeta, filesList]]                 // Output the single meta and list of files
        }
        .map { id, chr, group -> [chr, id, group] }
        .combine(ch_chroms.map { meta, bed -> [meta.chr, bed] }, by: 0) // add the chromosome bed files
        // move bed inside of the group
        .map { chr, id, group, bed -> [id, chr, [group[0], [bed, group[1]].flatten()]] }
        .set { ch_test1 }

    ch_bwaob
        .filter { meta, tab -> meta.control }
        .concat(ch_substracted) //s remove cpm
        .concat(ch_part_smooth.filter { meta, tab -> meta.minusinput }) // keep only minusinput
        .map { meta, tab -> [ meta.id, meta.chr, [meta, tab] ] }
        .groupTuple( by: [0, 1] )
        .map { id, chr, group ->
            def sortedGroup = group
                .sort { a, b ->
                    // Sorting criteria:
                    a[0].RFD <=> b[0].RFD ?:           // Sort by RFD (true vs false)
                    a[0].cpm <=> b[0].cpm ?:          // Sort by cpm (true vs false)
                    a[0].strand <=> b[0].strand    // Sort by strand (forward vs reverse)
                }

            def primaryMeta = sortedGroup[0][0] + ['minusinput':true]      // Take the first meta (after sorting)
            def filesList = sortedGroup*.get(1)      // Extract list of files

            [id, chr, [primaryMeta, filesList]]                 // Output the single meta and list of files
        }
        .map { id, chr, group -> [chr, id, group] }
        .combine(ch_chroms.map { meta, bed -> [meta.chr, bed] }, by: 0) // add the chromosome bed files
        // move bed inside of the group
        .map { chr, id, group, bed -> [id, chr, [group[0], [bed, group[1]].flatten()]] }
        .set { ch_test2 }


    ch_test1
        .concat(ch_test2)
        .map { chr, id, group -> [group[0], group[1]] }
        .set { ch_test_all }



    COLLECT_PARTITIONS_BY_CHROMOSOME (
        ch_test_all
    )
    ch_versions = ch_versions.mix(COLLECT_PARTITIONS_BY_CHROMOSOME.out.versions.first())


    // TODO: REMOVE print channel to file for debugging
    COLLECT_PARTITIONS_BY_CHROMOSOME.out.txt
        .map {
            meta, txt ->
                "${meta}\t${txt}"
        }
        .collectFile( name: 'ch_collect_part_by_chrom_out.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    ch_s_uniq_meta
        .concat(ch_i_uniq_meta)
        .unique { it -> it[0] }
        // remove chr from meta information
        .map {
            id, meta ->
                def meta_clone = meta.clone()
                meta_clone.remove('chr')
                [ id, meta_clone ]
        }
        .set { ch_si_uniq_meta_mod }

    // collect all chromosomes by meta - meta.chr
    COLLECT_PARTITIONS_BY_CHROMOSOME.out.txt
        .collectFile(newLine: true, sort: false) { //, storeDir: "${params.outdir}/${params.aligner}/mergedLibrary/scarseq/collect") {
            meta, txt ->
                // to do it per id and minusinput:
                [ "${meta.id}${meta.minusinput ? '.minusinput' : ''}.final.txt", txt ]
        }
        // re-add meta information based on meta.id in the filename
        .map {
            txt ->
                def id = txt.name.split("\\.")[0]
                [ id, txt ]
        }
        .combine(ch_si_uniq_meta_mod, by: 0)
        // flip meta and tab to have meta first
        .map {
            id, txt, meta ->
                [ meta, txt ]
        }
        .map {
            meta, txt ->
                def minusinput = txt.name.contains(".minusinput")
                def meta_clone = meta + ['minusinput':minusinput ]
                [ meta_clone, txt ]
        }
        .set { ch_collected }

    //
    // MODULE: Generate the final partition bedgraph
    //
    FINAL_PARTITION_BEDGRAPH (
        ch_collected
    )
    ch_versions = ch_versions.mix(FINAL_PARTITION_BEDGRAPH.out.versions.first())

    //
    // MODULE: Sort the final partition bedgraph
    //
    FILE_SORT_PARTITIONS (
        FINAL_PARTITION_BEDGRAPH.out.tmp,
        'bdg'
    )
    ch_versions = ch_versions.mix(FILE_SORT_PARTITIONS.out.versions.first())

    //
    // MODULE: Convert the final partition bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_PARTITIONS (
        FILE_SORT_PARTITIONS.out.sorted,
        ch_chrom_sizes.map { it[1] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_PARTITIONS.out.versions.first())


    FINAL_PARTITION_BEDGRAPH.out.txt
        .map {
            meta, txt ->
                // remove "_minusinput" from the id
                meta.minusinput ? [ meta.id.replaceAll("_minusinput", ""), meta, txt ] :
                meta.control ? [ meta.control, meta, txt ] :
                [ meta.id, meta, txt ]
        }
        .branch { id, meta, txt ->
            control:    !meta.control && !meta.minusinput
            samples:    meta.control && !meta.minusinput
            minusinput: meta.minusinput
        }
        .set { ch_part_to_plot }
    

    ch_part_to_plot.samples
        .combine(ch_part_to_plot.control, by: 0)
        .map {
            id, meta1, txt1, meta2, txt2 ->
                [ meta1.id, meta1, txt1, txt2 ]
        }
        .combine(ch_part_to_plot.minusinput, by: 0)
        .map {
            id, meta1, txt1, txt2, meta3, txt3 ->
                def okseq = meta1.okseq_part_file ? file(meta1.okseq_part_file) : null
                [ meta1, txt1, txt2, txt3, okseq ]
        }
        .set { ch_part_to_plot }

    //
    // MODULE: Plot the final partition
    //
    FINAL_PARTITION_PLOT (
        ch_part_to_plot,
        ch_blacklist,
        ch_initiation_zones
        //ch_scaffolds
    )
    ch_versions = ch_versions.mix(FINAL_PARTITION_PLOT.out.versions.first())


    emit:
    tab      = PARTITION_SMOOTH.out.rfd       // channel: [ val(meta), [ tab ] ]
    versions = ch_versions                    // channel: [ versions.yml ]
}

