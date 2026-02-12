process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c36472269e8898f63b7b65dd40433462d541f9e75f9401f0bf8488021275d006/data' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(splicesites)
    tuple val(meta4), path(fasta)
    val   save_unaligned
    val   sort_bam

    output:
    tuple val(meta), path("*.sam")      , emit: sam     , optional:true
    tuple val(meta), path("*.bam")      , emit: bam     , optional:true
    tuple val(meta), path("*.cram")     , emit: cram    , optional:true
    tuple val(meta), path("*.csi")      , emit: csi     , optional:true
    tuple val(meta), path("*.crai")     , emit: crai    , optional:true
    tuple val(meta), path("*.log")      , emit: summary
    tuple val(meta), path("*fastq.gz")  , emit: fastq   , optional:true
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n '1s/.*version //p'"), emit: versions_hisat2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    def unaligned = ""
    def reads_args = ""
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-U ${reads}"
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
    }

    def ss = splicesites ? "--known-splicesite-infile ${splicesites}" : ""
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension == "cram" ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    """
    mkdir ./tmp && chmod 777 ./tmp

    INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
    hisat2 \\
        -x \$INDEX \\
        ${reads_args} \\
        ${ss} \\
        --summary-file ${prefix}.hisat2.summary.log \\
        --threads ${task.cpus} \\
        ${unaligned} \\
        ${args} \\
        --temp-directory ./tmp \\
        | samtools ${samtools_command} ${args2} --threads ${task.cpus} ${reference} -o ${prefix}.${extension} -
   
    if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
        mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
    fi

    if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
        mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
    fi
    """
    
    stub:
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension = (args2 ==~ extension_pattern) ? (args2 =~ extension_pattern)[0][2].toLowerCase() : "bam"
    def create_unmapped = ""
    if (meta.single_end) {
        create_unmapped = save_unaligned ? "touch ${prefix}.unmapped.fastq.gz" : ""
    } else {
        create_unmapped = save_unaligned ? "touch ${prefix}.unmapped_1.fastq.gz && touch ${prefix}.unmapped_2.fastq.gz" : ""
    }
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    } else if (extension == "bam") {
        create_index = "touch ${prefix}.csi"
    }

    """
    touch ${prefix}.${extension}
    ${create_index}
    touch ${prefix}.hisat2.summary.log
    ${create_unmapped}
    """

}
