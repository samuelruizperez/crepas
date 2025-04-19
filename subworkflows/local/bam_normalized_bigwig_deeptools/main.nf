include { DEEPTOOLS_BAMCOVERAGE } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { FILE_SORT } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG     } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BAM_NORMALIZED_BIGWIG_DEEPTOOLS {

    take:
    ch_bam_bai               // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_chrom_sizes          // channel: [ bed ]
    genome                  // string: genome name
    spikein_genome          // string: spike-in genome name
    skip_srpm           // boolean: skip the SRPM normalization step
    skip_cisrpm          // boolean: skip the CISRPM normalization step
    skip_cisrpmsoi      // boolean: skip the CISRPM-SOI normalization step

    main:

    ch_versions = Channel.empty()


    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1
            meta_clone.norm_factor_type = 'raw'
            [ meta_clone, bam, bai ]
        }
        .set { ch_bam_bai_raw }


    ch_bam_bai_rpm = Channel.empty()

    ch_bam_bai
        .map {
            meta, bam, bai ->
                def meta_clone = meta.clone()
                meta_clone.norm_factor_val = 1e6 / meta_clone.total_mapped_reads
                meta_clone.norm_factor_type = 'rpm'
                [ meta_clone, bam, bai ]
        }
        .set { ch_bam_bai_rpm }


    ch_bam_bai_srpm = Channel.empty()
    if (!skip_srpm) {

        ch_bam_bai
            .map {
                meta, bam, bai ->
                    [ meta.id, meta, bam, bai ]
            }
            .branch { id, meta, bam, bai ->
                endo: meta.genome == genome
                exo: meta.genome == spikein_genome
            }
            .set { ch_bam_bai_genome }

        ch_bam_bai_genome
            .endo
            .combine(ch_bam_bai_genome.exo, by: 0)
            .map { id, endo_meta, endo_bam, endo_bai, exo_meta, exo_bam, exo_bai ->
                def meta_clone = endo_meta.clone()
                meta_clone.norm_factor_val = 1e6 / exo_meta.total_mapped_reads
                meta_clone.norm_factor_type = 'srpm'
                [ meta_clone, endo_bam, endo_bai ]
            }
            .set { ch_bam_bai_srpm }
    }

    ch_bam_bai_cisrpm = Channel.empty()
    if (!skip_cisrpm) {

        // Split BAMs by genome (endo and exo) and by type (ip and control)
        ch_bam_bai
            .map {
                meta, bam, bai ->
                    [ meta.id, meta, bam, bai ]
            }
            .branch { id, meta, bam, bai ->
                endo_ip: meta.genome == genome && !meta.is_control
                endo_control: meta.genome == genome && meta.is_control
                exo_ip: meta.genome == spikein_genome && !meta.is_control
                exo_control: meta.genome == spikein_genome && meta.is_control
            }
            .set { ch_bam_bai_genome_type }

        // Combine the endo and exo BAMs (ChIPs)
        ch_bam_bai
            .endo_ip
            .combine(ch_bam_bai_genome_type.exo_ip, by: 0)
            .map { 
                ip_id, endo_ip_meta, endo_ip_bam, endo_ip_bai, exo_ip_meta, exo_ip_bam, exo_ip_bai ->
                    [ endo_ip_meta.control, endo_ip_meta, endo_ip_bam, endo_ip_bai, exo_ip_meta, exo_ip_bam, exo_ip_bai ]
            }
            .set { ch_bam_bai_genome_ip }

        // Combine the endo and exo BAMs (inputs)
        ch_bam_bai
            .endo_control
            .combine(ch_bam_bai_genome_type.exo_control, by: 0)
            .map { control_id, endo_control_meta, endo_control_bam, endo_control_bai, exo_control_meta, exo_control_bam, exo_control_bai ->
                [ endo_control_meta.id, endo_control_meta, endo_control_bam, endo_control_bai, exo_control_meta, exo_control_bam, exo_control_bai ]
            }
            .set { ch_bam_bai_genome_control }

        // Combine the combined ChIPs with the combined inputs
        ch_bam_bai_genome_ip
            .combine(ch_bam_bai_genome_control, by: 0)
            .map { 
                id, endo_ip_meta, endo_ip_bam, endo_ip_bai, exo_ip_meta, exo_ip_bam, exo_ip_bai, endo_control_meta, endo_control_bam, endo_control_bai, exo_control_meta, exo_control_bam, exo_control_bai ->
                    def meta_clone = endo_ip_meta.clone()
                    meta_clone.norm_factor_val = (1e6 / exo_ip_meta.total_mapped_reads) * (exo_control_meta.total_mapped_reads / endo_control_meta.total_mapped_reads)
                    meta_clone.norm_factor_type = 'cisrpm'
                    [ meta_clone, endo_ip_bam, endo_ip_bai ]
            }
            .set { ch_bam_bai_ip_cisrpm }
        
        // Now do the missing CISRPM for the endogenous inputs
        // In this case CISRPM is the same as RPM
        ch_bam_bai_rpm
            .filter { meta, bam, bai ->
                meta.genome == genome && meta.is_control
            }
            .map { meta, bam, bai ->
                def meta_clone = meta.clone()
                // just change norm_factor_type to cisrpm but keeping norm_factor_val the untouched
                meta_clone.norm_factor_type = 'cisrpm'
                [ meta_clone, bam, bai ]
            }
            .set { ch_bam_bai_control_cisrpm }

        ch_bam_bai_cisrpm = ch_bam_bai_ip_cisrpm.mix(ch_bam_bai_control_cisrpm)

    }

    ch_bam_bai_raw
        .mix(ch_bam_bai_rpm)
        .mix(ch_bam_bai_srpm)
        .mix(ch_bam_bai_cisrpm)
        .set { ch_bam_bai_all }

    // TODO: print for debugging
    ch_bam_bai_all
        .map {
            meta, bam, bai ->
                "${meta}\t${bam}\t${bai}"
        }
        .collectFile( name: 'ch_bam_bai_all_norms.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    //
    // MODULE: Calculate raw and normalized coverage per bin
    //
    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai_all,
        [],
        []
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())

    //
    // MODULE: Sort the bedgraph
    //
    FILE_SORT (
       DEEPTOOLS_BAMCOVERAGE.out.bedgraph,
        'bedgraph'
    )
    ch_versions = ch_versions.mix(FILE_SORT.out.versions.first())

    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG (
        FILE_SORT.out.sorted,
        ch_chrom_sizes.map { it[1] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())
   

    emit:
    bigwig      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig        // channel: [ val(meta), [ bigwig ] ]

    versions      = ch_versions                            // channel: [ versions.yml ]
}

