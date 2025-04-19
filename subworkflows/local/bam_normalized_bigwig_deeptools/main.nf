
include { BAM_FLAGSTAT_MAPPED   } from '../../../modules/local/bam_flagstat_mapped/main'
include { DEEPTOOLS_BAMCOVERAGE } from '../../../modules/nf-core/deeptools/bamcoverage/main'


include { FILE_SORT } from '../../../modules/local/file_sort/main'

include { UCSC_BEDGRAPHTOBIGWIG     } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BAM_NORMALIZED_BIGWIG_DEEPTOOLS {

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

