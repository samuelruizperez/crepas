#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grothlab/glseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/grothlab/glseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GLSEQ                   } from './workflows/glseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_grothlab_glseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_grothlab_glseq_pipeline'
include { getMacsGsize            } from './subworkflows/local/utils_grothlab_glseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_grothlab_glseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GROTHLAB_GLSEQ {

    main:
    ch_versions = Channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        GENOME PARAMETER VALUES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // collect paths from genome attributes file (e.g. iGenomes.config; optional)
    // we cannot overwrite params in the workflow (they stay null as coming from the config file)
    // TODO: simplify, readability
    def fasta               = params.containsKey('fasta') ? params.fasta : (params.refgenie_ignore ? null : getGenomeAttribute('fasta'))
    def bwa_index           = params.containsKey('bwa_index') ? params.bwa_index : (params.refgenie_ignore ? null : getGenomeAttribute('bwa'))
    def bowtie2_index       = params.containsKey('bowtie2_index') ? params.bowtie2_index : (params.refgenie_ignore ? null : getGenomeAttribute('bowtie2_index'))
    def chromap_index       = params.containsKey('chromap_index') ? params.chromap_index : (params.refgenie_ignore ? null : getGenomeAttribute('chromap'))
    def star_index          = params.containsKey('star_index') ? params.star_index : (params.refgenie_ignore ? null : getGenomeAttribute('star'))
    def hisat2_index        = params.containsKey('hisat2_index') ? params.hisat2_index : (params.refgenie_ignore ? null : getGenomeAttribute('hisat2'))
    def gtf                 = params.containsKey('gtf') ? params.gtf : (params.refgenie_ignore ? null : getGenomeAttribute('gtf'))
    def gff                 = params.containsKey('gff') ? params.gff : (params.refgenie_ignore ? null : getGenomeAttribute('gff'))
    def gene_bed            = params.containsKey('gene_bed') ? params.gene_bed : (params.refgenie_ignore ? null : getGenomeAttribute('gene_bed'))
    def blacklist           = params.containsKey('blacklist') ? params.blacklist : (params.refgenie_ignore ? null : getGenomeAttribute('blacklist'))
    def sparsebed           = params.containsKey('sparsebed') ? params.sparsebed : (params.refgenie_ignore ? null : getGenomeAttribute('sparsebed'))
    def active_regions      = params.containsKey('active_regions') ? params.active_regions : (params.refgenie_ignore ? null : getGenomeAttribute('active_regions'))
    def rocco_params        = params.containsKey('rocco_params') ? params.rocco_params : (params.refgenie_ignore ? null : getGenomeAttribute('rocco_params'))
    def splicesites         = params.containsKey('splicesites') ? params.splicesites : (params.refgenie_ignore ? null : getGenomeAttribute('splicesites'))
    def initiation_zones    = params.containsKey('initiation_zones') ? params.initiation_zones : (params.refgenie_ignore ? null : getGenomeAttribute('initiation_zones'))
    def tecount_gene_index  = params.containsKey('tecount_gene_index') ? params.tecount_gene_index : (params.refgenie_ignore ? null : getGenomeAttribute('tecount_gene_index'))
    def telocal_gene_index  = params.containsKey('telocal_gene_index') ? params.telocal_gene_index : (params.refgenie_ignore ? null : getGenomeAttribute('telocal_gene_index'))
    def te_gtf              = params.containsKey('te_gtf') ? params.te_gtf : (params.refgenie_ignore ? null : getGenomeAttribute('te_gtf'))
    def tecount_te_index    = params.containsKey('tecount_te_index') ? params.tecount_te_index : (params.refgenie_ignore ? null : getGenomeAttribute('tecount_te_index'))
    def telocal_te_index    = params.containsKey('telocal_te_index') ? params.telocal_te_index : (params.refgenie_ignore ? null : getGenomeAttribute('telocal_te_index'))
    def macs_gsize          = params.containsKey('macs_gsize') ? params.macs_gsize : getMacsGsize(params)


    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME (
        params.genome,
        params.spikein_genome,
        params.aligner,
        fasta,
        gtf,
        gff,
        blacklist,
        params.read_length,
        macs_gsize,
        sparsebed,
        active_regions,
        rocco_params,
        gene_bed,
        bwa_index,
        bowtie2_index,
        chromap_index,
        star_index,
        hisat2_index,
        splicesites,
        initiation_zones,
        params.skip_te_counting,
        params.skip_telocal,
        tecount_gene_index,
        telocal_gene_index,
        te_gtf,
        tecount_te_index,
        telocal_te_index
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // WORKFLOW: Run grothlab/glseq workflow
    //
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))
    GLSEQ (
        ch_samplesheet,
        ch_versions,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.chrom_sizes_endo,
        PREPARE_GENOME.out.chrom_sizes_exo,
        PREPARE_GENOME.out.effective_gsize,
        PREPARE_GENOME.out.effective_gfraction,
        PREPARE_GENOME.out.whitelist,
        PREPARE_GENOME.out.blacklist,
        PREPARE_GENOME.out.sparsebed,
        PREPARE_GENOME.out.active_regions,
        PREPARE_GENOME.out.rocco_params,
        PREPARE_GENOME.out.initiation_zones,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.chromap_index,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.splicesites,
        PREPARE_GENOME.out.tecount_gene_index,
        PREPARE_GENOME.out.telocal_gene_index,
        PREPARE_GENOME.out.tecount_te_index,
        PREPARE_GENOME.out.telocal_te_index
    )

    emit:
    multiqc_report = GLSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [version1, version2, ...]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    GROTHLAB_GLSEQ ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        GROTHLAB_GLSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
