process BAM_SPLIT_BY_CHROMOSOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam), path(bai), val(chrom_list)
    val index_format

    output:
    tuple val(meta), path("${prefix}.bam"),                 emit: bam,  optional: true
    tuple val(meta), path("${prefix}.cram"),                emit: cram, optional: true
    tuple val(meta), path("${prefix}.sam"),                 emit: sam,  optional: true
    tuple val(meta), path("${prefix}.${extension}.crai"),   emit: crai, optional: true
    tuple val(meta), path("${prefix}.${extension}.csi"),    emit: csi,  optional: true
    tuple val(meta), path("${prefix}.${extension}.bai"),    emit: bai,  optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.split_by_chrom"
    def chrom_list_arg = chrom_list ? "${chrom_list.join(' ')}" : ''
    extension = args.contains("--output-fmt sam") ? "sam" :
                args.contains("--output-fmt cram") ? "cram" :
                "bam"
    output_file = index_format ? "${prefix}.${extension}##idx##${prefix}.${extension}.${index_format} --write-index" : "${prefix}.${extension}"
    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error "Index format not one of bai, csi, crai."
        } else if (extension == "sam") {
            error "Indexing not compatible with SAM output"
        }
    }
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools view \\
        ${args} \\
        --threads ${task.cpus - 1} \\
        -b ${bam} \\
        ${chrom_list_arg} \\
        -o ${output_file}
    """

    stub:
    def args  = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.split_by_chrom"
    extension = args.contains("--output-fmt sam") ? "sam" :
                args.contains("--output-fmt cram") ? "cram" :
                "bam"
    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error "Index format not one of bai, csi, crai."
        } else if (extension == "sam") {
            error "Indexing not compatible with SAM output"
        }
    }
    index = index_format ? "touch ${prefix}.${extension}.${index_format}" : ""
    """
    touch  ${prefix}.${extension}
    ${index}
    """
}
