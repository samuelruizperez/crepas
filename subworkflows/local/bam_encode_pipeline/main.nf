include { SAMTOOLS_SORT               } from '../../../modules/nf-core/samtools/sort/main'
include { BEDTOOLS_BAMTOBED           } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { BED_TO_TAGALIGN            } from '../../../modules/local/bed_to_tagalign/main'
include { TAGALIGN_SELF_PSEUDOREPLICATES } from '../../../modules/local/tagalign_self_pseudoreplicates/main'
include { CAT_CAT as TAGALIGN_POOL } from '../../../modules/nf-core/cat/cat/main'
include { PHANTOMPEAKQUALTOOLS } from '../../../modules/nf-core/phantompeakqualtools/main'


workflow BAM_ENCODE_PIPELINE {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_fasta                          // channel: [ fasta ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Name-sorting BAM files
    //

    SAMTOOLS_SORT(
        ch_bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())


    // 2a section of the ENCODE 3 ChIP-seq pipeline:

    //
    // MODULE: Convert BAM to BED
    //
    BEDTOOLS_BAMTOBED (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())

    //
    // MODULE: Convert BED to TAGALIGN
    //
    BED_TO_TAGALIGN (
        BEDTOOLS_BAMTOBED.out.bed
    )
    ch_versions = ch_versions.mix(BED_TO_TAGALIGN.out.versions.first())

    //
    // MODULE: Generate self-pseudoreplicates from TAGALIGN files
    //
    TAGALIGN_SELF_PSEUDOREPLICATES (
        BED_TO_TAGALIGN.out.tagalign
    )
    ch_versions = ch_versions.mix(TAGALIGN_SELF_PSEUDOREPLICATES.out.versions.first())

    // Add pseudoreplicate to metadata
    TAGALIGN_SELF_PSEUDOREPLICATES.out.tagalign1
        .map { meta, tagalign -> [ meta + [ pseudoreplicate: '1' ], tagalign ]}
        .mix(
            TAGALIGN_SELF_PSEUDOREPLICATES.out.tagalign2
                .map { meta, tagalign -> [ meta + [ pseudoreplicate: '2' ], tagalign ] }
        )
        .set {ch_self_pseudoreps}
    

    BED_TO_TAGALIGN.out.tagalign
        .mix(ch_self_pseudoreps)
        .map { meta, tagalign ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            def antibody = meta.input_control_of_antibody ?: meta.antibody
            def pseudoreplicate = meta.pseudoreplicate ?: false
            [ meta.id, antibody, pseudoreplicate, meta, tagalign ]
        }
        .groupTuple(by: [0, 1, 2])
        .map { id, antibody, pseudoreplicate, metas, tagaligns ->
            [ metas[0], tagaligns.flatten() ]
        }
        .set { ch_tas_reps_and_pseudoreps }

    //
    // MODULE: Pool replicates and pseudoreplicates with cat
    //
    TAGALIGN_POOL (
        ch_tas_reps_and_pseudoreps
    )
    ch_versions = ch_versions.mix(TAGALIGN_POOL.out.versions.first())

    //
    // MODULE: Call peaks with phantompeakqualtools SPP
    //
    // PHANTOMPEAKQUALTOOLS (
    //     TAGALIGN_POOL.out.file_out
    // )
    // ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions.first())


    emit:


    versions                     = ch_versions                      // channel: [ versions.yml ]
}
