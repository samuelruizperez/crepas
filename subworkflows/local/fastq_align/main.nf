//
// FASTQ_ALIGN: Align fastq files to a reference genome
//

include { FASTQ_ALIGN_BWA as FASTQ_ALIGN_BWAMEM1     } from '../../../subworkflows/nf-core/fastq_align_bwa'
include { FASTQ_ALIGN_BWAMEM2                        } from '../../../subworkflows/local/fastq_align_bwamem2'
include { FASTQ_ALIGN_BOWTIE2                        } from '../../../subworkflows/nf-core/fastq_align_bowtie2'
include { FASTQ_ALIGN_BOWTIE                         } from '../../../subworkflows/local/fastq_align_bowtie'
include { FASTQ_ALIGN_STROBEALIGN                    } from '../../../subworkflows/local/fastq_align_strobealign'
include { FASTQ_ALIGN_CHROMAP                        } from '../../../subworkflows/nf-core/fastq_align_chromap'
include { FASTQ_ALIGN_STAR                           } from '../../../subworkflows/nf-core/fastq_align_star'
include { FASTQ_ALIGN_HISAT2                         } from '../../../subworkflows/nf-core/fastq_align_hisat2'
include { FASTQ_ALIGN_MINIMAP2                       } from '../../../subworkflows/local/fastq_align_minimap2'

workflow FASTQ_ALIGN {
    take:
    ch_reads            // channel: [mandatory] meta, reads
    ch_fasta
    aligners
    ch_bwa_index
    ch_bwamem2_index
    ch_bowtie_index
    ch_bowtie2_index
    ch_chromap_index
    ch_star_index
    ch_hisat2_index
    ch_minimap2_index
    ch_gtf
    ch_splicesites
    save_unaligned
    seq_platform
    seq_center
    sort_bam

    main:

    ch_multiqc_files          = channel.empty()
    ch_samtools_stats_summary = channel.empty()  
    ch_versions               = channel.empty()
    ch_bam                    = channel.empty()
    ch_bai                    = channel.empty()

    if (aligners.contains('bwa')) {
        FASTQ_ALIGN_BWAMEM1 (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'bwa' ], reads ] },
            ch_bwa_index,
            false,
            ch_fasta
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_BWAMEM1.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_BWAMEM1.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BWAMEM1.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM1.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM1.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM1.out.idxstats.collect { it -> it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BWAMEM1.out.versions.first())
    }

    if (aligners.contains('bwamem2')) {
        FASTQ_ALIGN_BWAMEM2 (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'bwamem2' ], reads ] },
            ch_bwamem2_index,
            false,
            ch_fasta
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_BWAMEM2.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_BWAMEM2.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BWAMEM2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM2.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM2.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM2.out.idxstats.collect { it -> it[1] })
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie & BAM QC
    //
    if (aligners.contains('bowtie')) {
        FASTQ_ALIGN_BOWTIE (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'bowtie' ], reads ] },
            ch_bowtie_index,
            save_unaligned,
            ch_fasta
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_BOWTIE.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_BOWTIE.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BOWTIE.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.idxstats.collect { it -> it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie2 & BAM QC
    //
    if (aligners.contains('bowtie2')) {
        FASTQ_ALIGN_BOWTIE2 (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'bowtie2' ], reads ] },
            ch_bowtie2_index,
            save_unaligned,
            sort_bam,
            ch_fasta
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_BOWTIE2.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_BOWTIE2.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BOWTIE2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.idxstats.collect { it -> it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with strobealign
    //
    if (aligners.contains('strobealign')) {
        FASTQ_ALIGN_STROBEALIGN (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'strobealign' ], reads ] },
            channel.value([[:], []]),
            ch_fasta,
            true
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_STROBEALIGN.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_STROBEALIGN.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_STROBEALIGN.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STROBEALIGN.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STROBEALIGN.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STROBEALIGN.out.idxstats.collect { it -> it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_STROBEALIGN.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with Chromap & BAM QC
    //
    if (aligners.contains('chromap')) {
        FASTQ_ALIGN_CHROMAP (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'chromap' ], reads ] },
            ch_chromap_index,
            ch_fasta,
            [],
            [],
            [],
            []
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_CHROMAP.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_CHROMAP.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_CHROMAP.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.idxstats.collect { it -> it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_CHROMAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    if (aligners.contains('star')) {
        FASTQ_ALIGN_STAR (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'star' ], reads ] },
            ch_star_index,
            ch_gtf,
            true,
            seq_platform ?: '',
            seq_center ?: '',
            ch_fasta,
            channel.value([[:], []])
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_STAR.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_STAR.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_STAR.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.idxstats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.log_final.collect { it -> it[1] })
    }

    if (aligners.contains('hisat2')) {
        FASTQ_ALIGN_HISAT2 (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'hisat2' ], reads ] },
            ch_hisat2_index,
            ch_splicesites,
            ch_fasta,
            save_unaligned,
            sort_bam
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_HISAT2.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_HISAT2.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_HISAT2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary.collect { it -> it[1] })
    }

    if (aligners.contains('minimap2')) {
        FASTQ_ALIGN_MINIMAP2 (
            ch_reads.map { meta, reads -> [ meta + [ aligner: 'minimap2' ], reads ] },
            ch_minimap2_index,
            ch_fasta,
            true,
            'bai',
            false,
            true
        )
        ch_bam = ch_bam.mix(FASTQ_ALIGN_MINIMAP2.out.bam)
        ch_bai = ch_bai.mix(FASTQ_ALIGN_MINIMAP2.out.bai)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_MINIMAP2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIMAP2.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIMAP2.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIMAP2.out.idxstats.collect { it -> it[1] })
    }

    // Prepend aligner name to meta.id
    ch_bam
        .map { meta, bam ->
            def meta_clone = meta.clone()
            meta_clone.id = "${meta.aligner}_${meta.id}"
            [ meta_clone, bam ]
        }
        .set { ch_bam }
    
    ch_bai
        .map { meta, bai ->
            def meta_clone = meta.clone()
            meta_clone.id = "${meta.aligner}_${meta.id}"
            [ meta_clone, bai ]
        }
        .set { ch_bai }

    emit:
    bam      = ch_bam                                   // channel: [ val(meta), path(bam) ]
    bai      = ch_bai                                   // channel: [ val(meta), path(bai) ]
    samtools_stats_summary = ch_samtools_stats_summary  // channel: [ val(meta), path(summary) ]
    multiqc_files = ch_multiqc_files                    // channel: [ path(multiqc files ]
    versions = ch_versions                              // channel: [ path(versions.yml) ]
}