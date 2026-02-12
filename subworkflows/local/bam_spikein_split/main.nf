//
// Splitting input BAMs by spike-in genome (into endogenous and exogenous BAMs)
//

include { BAM_SPLIT_BY_GENOME as BAM_SPLIT_BY_GENOME_ENDO } from '../../../modules/local/bam_split_by_genome/main'
include { BAM_SPLIT_BY_GENOME as BAM_SPLIT_BY_GENOME_EXO } from '../../../modules/local/bam_split_by_genome/main'
include { SAMTOOLS_SORT           } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_REMOVE_ORPHANS      } from '../../../modules/local/bam_remove_orphans/main'
include { BAM_SORT_STATS_SAMTOOLS  } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow BAM_SPIKEIN_SPLIT {
    take:
    ch_bam               // channel: [ val(meta), [ bam ], [bai] ]
    ch_fasta             // channel: [ val(meta), path(fasta) ]
    genome               // val
    spikein_genome       // val

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: split BAMs by spike-in genome (keep endogenous)
    //
    BAM_SPLIT_BY_GENOME_ENDO (
        ch_bam,
        genome,
        spikein_genome,
        'endo'
    )
    ch_versions = ch_versions.mix(BAM_SPLIT_BY_GENOME_ENDO.out.versions.first())

    //
    // MODULE: split BAMs by spike-in genome (keep exogenous)
    //
    BAM_SPLIT_BY_GENOME_EXO (
        ch_bam,
        genome,
        spikein_genome,
        'exo'
    )
    ch_versions = ch_versions.mix(BAM_SPLIT_BY_GENOME_EXO.out.versions.first())

    // add genome as meta field
    ch_bam_endo = BAM_SPLIT_BY_GENOME_ENDO.out.bam.map { it -> [ it[0] + [ genome: genome ], it[1] ] }
    ch_bam_exo = BAM_SPLIT_BY_GENOME_EXO.out.bam.map { it -> [ it[0] + [ genome: spikein_genome ], it[1] ] }
    ch_bam = ch_bam_endo.mix(ch_bam_exo)

    // Separate single-end and paired-end BAM files (SE do not have orphans)
    ch_bam
        .branch { meta, bam ->
            se: meta.single_end
            pe: !meta.single_end
        }
        .set { ch_bam }

    //
    // MODULE: Sort BAM files by query name
    //
    SAMTOOLS_SORT (
        ch_bam.pe,
        ch_fasta,
        ''
    )
    ch_bam_pe = SAMTOOLS_SORT.out.bam

    //
    // MODULE: Remove orphan reads left by splitting by genome (only for paired-end BAMs)
    //
    BAM_REMOVE_ORPHANS (
        ch_bam_pe,
        true
    )
    ch_versions = ch_versions.mix(BAM_REMOVE_ORPHANS.out.versions.first())

    // Mixing SE and PE (removed orphans) files
    ch_filtered_bam = ch_bam.se.mix(BAM_REMOVE_ORPHANS.out.bam)

    //
    // MODULE: Sort BAM by coordinate and generate stats
    //
    // TODO: I would need to separate the genome fasta, right now both endo and exo stats
    // are analyzed with the same (main) genome
    // either way, this has no effect on BAM output (https://bioinformatics.stackexchange.com/a/4218)
    BAM_SORT_STATS_SAMTOOLS (
        ch_filtered_bam,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.stats.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it -> it[1] })

    emit:

    bam           = BAM_SORT_STATS_SAMTOOLS.out.bam
    endo_bam      = BAM_SORT_STATS_SAMTOOLS.out.bam.filter { it -> it[0].genome == genome }            // channel: [ val(meta), [ bam ] ]
    exo_bam       = BAM_SORT_STATS_SAMTOOLS.out.bam.filter { it -> it[0].genome == spikein_genome }    // channel: [ val(meta), [ bam ] ]

    bai           = BAM_SORT_STATS_SAMTOOLS.out.bai
    endo_bai      = BAM_SORT_STATS_SAMTOOLS.out.bai.filter { it -> it[0].genome == genome }            // channel: [ val(meta), [ bai ] ]
    exo_bai       = BAM_SORT_STATS_SAMTOOLS.out.bai.filter { it -> it[0].genome == spikein_genome }

    stats         = BAM_SORT_STATS_SAMTOOLS.out.stats
    endo_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats.filter { it -> it[0].genome == genome }            // channel: [ val(meta), [ stats ] ]
    exo_stats     = BAM_SORT_STATS_SAMTOOLS.out.stats.filter { it -> it[0].genome == spikein_genome }    // channel: [ val(meta), [ stats ] ]

    flagstat      = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    endo_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat.filter { it -> it[0].genome == genome }            // channel: [ val(meta), [ flagstat ] ]
    exo_flagstat  = BAM_SORT_STATS_SAMTOOLS.out.flagstat.filter { it -> it[0].genome == spikein_genome }    // channel: [ val(meta), [ flagstat ] ]

    idxstats      = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    endo_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats.filter { it -> it[0].genome == genome }            // channel: [ val(meta), [ idxstats ] ]
    exo_idxstats  = BAM_SORT_STATS_SAMTOOLS.out.idxstats.filter { it -> it[0].genome == spikein_genome }    // channel: [ val(meta), [ idxstats ] ]

    multiqc_files = ch_multiqc_files                      // channel: [ multiqc_files ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
