include { BEDTOOLS_MAKEWINDOWS           } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BED_SPLIT_BY_CHROMOSOME       } from '../../../modules/local/bed_split_by_chromosome/main'
include { BED_SPLIT_WINDOWS     } from '../../../modules/local/bed_split_windows/main'
include { UCSC_BIGWIGAVERAGEOVERBED  } from '../../../modules/nf-core/ucsc/bigwigaverageoverbed/main'
include { CPM_CALCULATION as CPM_CALCULATION_SAMPLES  } from '../../../modules/local/cpm_calculation/main'
include { CPM_CALCULATION as CPM_CALCULATION_INPUTS  } from '../../../modules/local/cpm_calculation/main'
include { NORMALIZE_STRANDS     } from '../../../modules/local/normalize_strands/main'

workflow SCAR_SMOOTH_PARTITIONS {

    take:
    ch_bigwig          // channel: [ val(meta), [ bigwig ] ]
    ch_chrom_sizes  // channel: [ bed ]

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

    // create channel by combining the bigwig files with the chromosomes (all chroms per bigwig)
    ch_bigwig
        .combine( BED_SPLIT_BY_CHROMOSOME.out.bed )
        .set { ch_bigwig_chroms }

    // here we have to separate the channel again because UCSC_BIGWIGAVERAGEOVERBED
    // expects a channel with bigwig files and a channel with chromosome beds
    ch_bw_combs = ch_bigwig_chroms.map {
        bigwig, chrom ->
            [ bigwig ]
    }

    // add bigwig meta information to the chromosome bed files
    ch_chroms_combs = ch_bigwig_chroms.map {
        bigwig, chrom ->
            [ bigwig.meta, chrom.map{ it[1] } ]
        }

    // Aggregate counts in windows
    UCSC_BIGWIGAVERAGEOVERBED (
        ch_chroms_combs,
        ch_bw_combs.map{ it[1] }
    )
    ch_bwaob = UCSC_BIGWIGAVERAGEOVERBED.out.tab
    ch_versions = ch_versions.mix(UCSC_BIGWIGAVERAGEOVERBED.out.versions.first())

    // separate tab files in samples and input controls


    // separate samples and input controls ( meta.control contains something or not)
    ch_bwaob
        .map {
            meta, chroms ->
                meta.control ? null : [ meta.id, chroms ]
        }
        .set { ch_inputs }

    ch_bwaob
        .map {
            meta, chroms ->
                meta.control ? [ meta.id, chroms ] : null
        }
        .set { ch_samples }


    // concatenating the sample and input tabs

    ch_samples
        .map { it[1] }
        .collectFile( name: 'ch_samples.bed', newLine: true, sort: false, storeDir: "${params.outdir}" )
        .set { ch_samples_bed }

    ch_inputs
        .map { it[1] }
        .collectFile( name: 'ch_inputs.bed', newLine: true, sort: false, storeDir: "${params.outdir}" )
        .set { ch_inputs_bed }

    // CPM calculation
    CPM_CALCULATION_SAMPLES (
        ch_samples_bed
    )
    ch_versions = ch_versions.mix(CPM_CALCULATION_SAMPLES.out.versions.first())

    CPM_CALCULATION_INPUTS (
        ch_inputs_bed
    )
    ch_versions = ch_versions.mix(CPM_CALCULATION_INPUTS.out.versions.first())


    // concat samples and inputs and add their corresponding cpm
    ch_samples
        .map {
            meta, tab ->
                [ meta, tab, CPM_CALCULATION_SAMPLES.out.cpm[meta.id].splitCsv(header:false) ]
        }
        .concat(

            ch_inputs
            .map {
                meta, tab ->
                    [ meta, tab, CPM_CALCULATION_INPUTS.out.cpm[meta.id].splitCsv(header:false) ]
            }
        )
        .set { ch_samples_inputs_cpm }

    NORMALIZE_STRANDS (
        ch_samples_inputs_cpm,
    )


   emit:
   // bigwig      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig   // channel: [ val(meta), [ bigwig ] ]

      versions = ch_versions                         // channel: [ versions.yml ]
}

