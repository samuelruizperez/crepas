//
// SUBWORKFLOW: Calculating RRPM values from CHOR-seq BAM files as described in:
// https://www.nature.com/articles/s41596-021-00585-3#Sec60
//

include { CHOR_NORM_FACTOR_CALCULATION                                     } from '../../../modules/local/chor_norm_factor_calculation/main'
include { DEEPTOOLS_BAMCOVERAGE                                         } from '../../../modules/nf-core/deeptools/bamcoverage/main'
// include { BEDTOOLS_BAMTOBED                                         } from '../../../modules/nf-core/bedtools/bamtobed/main'
// include { BEDTOOLS_SLOP                                             } from '../../../modules/nf-core/bedtools/slop/main'
// include { BEDTOOLS_MAKEWINDOWS                                      } from '../../../modules/nf-core/bedtools/makewindows/main'
// include { BEDTOOLS_GENOMECOV                                       } from '../../../modules/nf-core/bedtools/genomecov/main'
include { FILE_SORT } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG     } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BAM_CHORSEQ_RRPM {

    take:
    ch_bam_bai               // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_chrom_sizes          // channel: [ bed ]
    genome                  // string: genome name
    spikein_genome          // string: spike-in genome name

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: count the total number of unique reads in the spike-in BAM file
    //
    ch_bam_bai
        .branch { meta, bam, bai ->
            endo: meta.genome == genome
            exo: meta.genome == spikein_genome
        }
        .set { ch_bam_bai_genome }
    
    CHOR_NORM_FACTOR_CALCULATION ( ch_bam_bai_genome.exo.map { it -> [ it[0], it[1] ] } )
    ch_versions = ch_versions.mix(CHOR_NORM_FACTOR_CALCULATION.out.versions.first())

    CHOR_NORM_FACTOR_CALCULATION.out.txt
        .map {
            meta, spikein_scale ->
                [ meta.id, spikein_scale.splitCsv(header:false)[0][0] ]
        }
        .set { ch_spikein_scale }


    // ch_spikein_scale
    //         .map {
    //             meta, scale ->
    //                 "${meta}\t${scale}"
    //         }
    //         .collectFile( name: 'ch_spikein_scale.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )
        
    // BEDTOOLS_BAMTOBED ( ch_bam_bai_genome.endo )
    // ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())

    // BEDTOOLS_SLOP (
    //     BEDTOOLS_BAMTOBED.out.bed,
    //     ch_chrom_sizes.map { it[1] }
    // )
    // ch_versions = ch_versions.mix(BEDTOOLS_SLOP.out.versions.first())

    // Add the spike-in normalization factor
    // Creating channel: [ val(meta), [ bed ], [ scale ] ] 
    ch_bam_bai_genome
        .endo
        .map {
            meta, bam, bai ->
                [ meta.id, meta, bam, bai ]
        }
        .combine(ch_spikein_scale, by: 0)
        .map {
            id, meta, bam, bai, scale ->
                // add scale to the meta
                def meta_clone = meta.clone()
                meta_clone.rrpm_scale = scale.toDouble()
                [ meta_clone, bam, bai ]
        }
        .set { ch_bam_bai_scale }

    ch_bam_bai_scale
                .map {
                    meta, bam, bai ->
                        "${meta}\t${bam}\t${bai}"
                }
                .collectFile( name: 'ch_bam_bai_scale.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )
            
    // BEDTOOLS_MAKEWINDOWS ( ch_chrom_sizes )
    // ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())

    // BEDTOOLS_GENOMECOV (
    //     ch_slop_scale,
    //     BEDTOOLS_MAKEWINDOWS.out.bed.map { it[1] },
    //     'bedGraph',
    //     true
    // )
    // ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai_scale,
        [],
        []
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())

    //
    // MODULE: Sort the final partition bedgraph
    //
    FILE_SORT (
       DEEPTOOLS_BAMCOVERAGE.out.bedgraph,
        'bedgraph'
    )
    ch_versions = ch_versions.mix(FILE_SORT.out.versions.first())


    UCSC_BEDGRAPHTOBIGWIG (
        FILE_SORT.out.sorted,
        ch_chrom_sizes.map { it[1] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())
   
    emit:
    rrpm      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig        // channel: [ val(meta), [ bigwig ] ]

    versions      = ch_versions                            // channel: [ versions.yml ]
}

