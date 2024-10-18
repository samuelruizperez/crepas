include { BEDTOOLS_MAKEWINDOWS           } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BED_SPLIT_WINDOWS     } from '../../../modules/local/bed_split_windows/main'
include { UCSC_BIGWIGAVERAGEOVERBED  } from '../../../modules/nf-core/ucsc/bigwigaverageoverbed/main'
include { NORMALIZE_STRANDS     } from '../../../modules/local/normalize_strands/main'

workflow SCAR_SMOOTH_PARTITIONS {

    take:
    //ch_bam          // channel: [ val(meta), [ bam ] ]
    //ch_chrom_sizes  // channel: [ bed ]

    main:

    ch_versions = Channel.empty()

    BEDTOOLS_MAKEWINDOWS (
            ch_chrom_sizes
        )
        ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())


        // TODO: for this step, spikein chromosomes should be excluded
        // creating a channel with each chromosome to iterate over
        ch_chrom_sizes
            .map {
                meta, bed ->
                    bed.splitCsv(header:false, sep:'\t')
            }
            .map{ it -> it[0] }
            .set { ch_chroms }

        // Split windows on chromosome
        BED_SPLIT_BY_CHROMOSOME (
            ch_chroms,
            BEDTOOLS_MAKEWINDOWS.out.bed.first()
        )
        ch_versions = ch_versions.mix(BED_SPLIT_BY_CHROMOSOME.out.versions.first())

        // combine the bigwig files with the chromosomes
        UCSC_BEDGRAPHTOBIGWIG.out.bigwig
            .combine( BED_SPLIT_BY_CHROMOSOME.out.bed )
            .set { ch_bigwig_chroms }

        ch_bw_combs = ch_bigwig_chroms.map {
            bigwig, chrom ->
                [ bigwig ]
        }

        ch_chroms_combs = ch_bigwig_chroms.map {
            bigwig, chrom ->
                [ chrom ]
        }

        // Aggregate counts in windows
        UCSC_BIGWIGAVERAGEOVERBED (
            ch_chroms_combs,
            ch_bw_combs.map{ it[1] }
        )
        ch_versions = ch_versions.mix(UCSC_BIGWIGAVERAGEOVERBED.out.versions.first())

        // separate tab files in samples and input controls


   emit:
   // bigwig      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig   // channel: [ val(meta), [ bigwig ] ]

      versions = ch_versions                         // channel: [ versions.yml ]
}

