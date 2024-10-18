include { BAM_SPLIT_BY_STRAND   } from '../../../modules/local/bam_split_by_strand/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BEDTOOLS_GENOMECOV             } from '../../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_SLOP                  } from '../../../modules/nf-core/bedtools/slop/main'
include { UCSC_BEDCLIP               } from '../../../modules/nf-core/ucsc/bedclip/main'
include { FILE_SORT                  } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG      } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'


workflow SCAR_CREATE_PARTITIONS {

    take:
    ch_bam          // channel: [ val(meta), [ bam ] ]
    ch_chrom_sizes  // channel: [ bed ]

    main:

    ch_versions = Channel.empty()

    BAM_SPLIT_BY_STRAND ( ch_bam )
    ch_versions = ch_versions.mix(BAM_SPLIT_BY_STRAND.out.versions.first())

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

    // Adding val(scale) to the input of the GENOMECOV module
    ch_bam
        .map {
            meta, bam ->
                [ meta, bam, 1 ]
        }
        .set { ch_bam_scale }

    // TODO: input channels for genomecov
    BEDTOOLS_GENOMECOV (
        ch_bam_scale,
        ch_chrom_sizes,
        'bdg',
        false
    )
    ch_genomecov = GENOMECOV.out.genomecov
    ch_versions  = ch_versions.mix(GENOMECOV.out.versions.first())

    BEDTOOLS_SLOP (
        ch_genomecov,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SLOP.out.versions.first())

    UCSC_BEDCLIP (
        BEDTOOLS_SLOP.out.bed,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions.first())

    // TODO: maybe a whole module for this is overkill
    SORT_FILE (
        UCSC_BEDCLIP.out.bed,
        'clip'
    )
    ch_versions = ch_versions.mix(SORT_FILE.out.versions.first())

    UCSC_BEDGRAPHTOBIGWIG (
        SORT_FILE.out.sorted,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    emit:
    bigwig      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig   // channel: [ val(meta), [ bigwig ] ]

      versions = ch_versions                         // channel: [ versions.yml ]
}

