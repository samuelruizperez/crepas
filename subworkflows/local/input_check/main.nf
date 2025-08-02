//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../../modules/local/samplesheet_check/main'
include { samplesheetToList                } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    seq_center  // string: sequencing center for read group

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it, seq_center) }
        .set { reads }

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(samplesheet, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2, fastq_umi ->
                def meta_clone = meta.clone()
                if (!fastq_2) {
                    meta_clone.single_end = true
                    return [ meta.id, meta.brep, meta_clone, [ fastq_1 ] ]
                } else {
                    meta_clone.single_end = false
                    return [ meta.id, meta.brep, meta_clone, [ fastq_1, fastq_2 ] ]
                }
                if (!fastq_umi) {
                    meta_clone.sep_umi_fq = false
                    return [ meta.id, meta.brep, meta_clone, [ fastq_1, fastq_2 ] ]
                } else {
                    meta_clone.sep_umi_fq = true
                    return [ meta.id, meta.brep, meta_clone, [ fastq_1, fastq_2, fastq_umi ] ]
                }
        }
        .set { ch_fastq }
    

    //
    // PARSING EXTRA META FIELDS
    //

    // get list of input control samples
    ch_fastq
        .map { it -> it[2].input_control }
        .unique()
        .toList()
        .map { it -> [it] } // unflatten
        .set { ch_ipcontrol_list }
    
    //ch_ipcontrol_list.view()


    // Add meta.is_input_control to ch_fastq
    ch_fastq
        //.combine(ch_ipcontrol_list)
        .map { id, brep, meta, fastqs -> //, ipcontrol_list ->
            def meta_clone = meta.clone()
            meta_clone.is_input_control = ch_ipcontrol_list.contains(meta.id)
            [ id, brep, meta_clone, fastqs ]
        }
        .set { ch_fastq }
        
    //ch_fastq.view()

    // Renaming ID with biological and technical replicate
    ch_fastq
        .map { id, brep, meta, fastqs ->
            def meta_clone = meta.clone()
            meta_clone.id = "${meta.exp_type}_${meta.id}_bRep_${brep}_tRep_${meta.trep}"

            def read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_tRep_.*$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\'"
            if (seq_center) {
                read_group = "\'@RG\\tID:${meta_clone.id}\\tSM:${meta_clone.id - ~/_tRep_.*$/}\\tPL:ILLUMINA\\tLB:${meta_clone.id}\\tPU:1\\tCN:${seq_center}\'"
            }
            meta_clone.read_group = read_group
            [ meta_clone.id, brep, meta_clone, fastqs ]
        }
        .set { ch_fastq }
        
    //ch_fastq.view()

    //
    // Checking consistency of meta fields after grouping
    //

    // get list of samples
    ch_fastq
        .map { it -> it[2].id }
        .unique()
        .toList()
        .map { it -> [it] } // unflatten
        .set { ch_sample_list }

    // Check that all input samples exist in meta.id
    ch_fastq
        .filter { it[2].is_input_control }
        .map { it -> it[2].id }



    // Group technical replicates
    ch_fastq
        .groupTuple(by: [0,1])
        .map { id, brep, metas, fastqs ->
            def 




        .map { samplesheet ->
            checkSamplesAfterGrouping(samplesheet)
        }
        .set { ch_fastq_grouped }


    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

//
// Function to check samples are internally consistent after being grouped
//
def checkSamplesAfterGrouping(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same strandedness
    def strandedness_ok = metas.collect{ it.strandedness }.unique().size == 1
    if (!strandedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must have the same strandedness!: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2, fastq_umi ] ]
def create_fastq_channel(LinkedHashMap row, String seq_center) {
    def meta = [:]
    meta.id                 = row.sample
    meta.single_end         = row.single_end.toBoolean()
    meta.sep_umi_fq         = row.sep_umi_fq.toBoolean()
    meta.okseq_part_file    = row.okseq_part_file
    meta.exp_type           = row.exp_type
    meta.strandedness       = row.strandedness
    meta.antibody           = row.antibody
    meta.control            = row.control
    meta.is_input_control         = row.is_input_control.toBoolean()

    def read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_T\d+$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\'"
    if (seq_center) {
        read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_T\d+$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\\tCN:${seq_center}\'"
    }
    meta.read_group = read_group

    // add path(s) of the fastq file(s) to the meta map
    // TODO: simplify this (fastq_umi)
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        if (meta.sep_umi_fq) {
            if (!file(row.fastq_umi).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> UMI FastQ file does not exist!\n${row.fastq_umi}"
            }
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_umi) ] ]
        } else {
            fastq_meta = [ meta, [ file(row.fastq_1) ] ]
        }
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        if (meta.sep_umi_fq) {
            if (!file(row.fastq_umi).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> UMI FastQ file does not exist!\n${row.fastq_umi}"
            }
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2), file(row.fastq_umi) ] ]
        } else {
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
        }
    }
    if (meta.okseq_part_file && !file(meta.okseq_part_file).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> OKSeq part file does not exist!\n${meta.okseq_part_file}"
    }
    return fastq_meta
}

