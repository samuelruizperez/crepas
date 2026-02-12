process HISAT2_BUILD {
    tag "${fasta}"
    label 'process_high'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c36472269e8898f63b7b65dd40433462d541f9e75f9401f0bf8488021275d006/data' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)
    tuple val(meta3), path(splicesites)

    output:
    tuple val(meta), path("hisat2"), emit: index
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n 's/.*version \\([^ ]*\\).*/\\1/p'"), emit: versions_hisat2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 0
    if (!task.memory) {
        log.info("[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this.")
    }
    else {
        log.info("[HISAT2 index build] Available memory: ${task.memory}")
        avail_mem = task.memory.toGiga()
    }

    def ss = ''
    def exon = ''
    def extract_exons = ''
    def hisat2_build_memory = params.hisat2_build_memory ? (params.hisat2_build_memory as MemoryUnit).toGiga() : 0
    if (avail_mem >= hisat2_build_memory) {
        log.info("[HISAT2 index build] At least ${hisat2_build_memory} GB available, so using splice sites and exons to build HISAT2 index")
        extract_exons = gtf ? "hisat2_extract_exons.py ${gtf} > ${gtf.baseName}.exons.txt" : ""
        ss = splicesites ? "--ss ${splicesites}" : ""
        exon = gtf ? "--exon ${gtf.baseName}.exons.txt" : ""
    }
    else {
        log.info("[HISAT2 index build] Less than ${hisat2_build_memory} GB available, so NOT using splice sites and exons to build HISAT2 index.")
        log.info("[HISAT2 index build] Use --hisat2_build_memory [small number] to skip this check.")
    }
    """
    mkdir hisat2
    ${extract_exons}
    hisat2-build \\
        -p ${task.cpus} \\
        ${ss} \\
        ${exon} \\
        ${args} \\
        ${fasta} \\
        hisat2/${fasta.baseName}
    """

    stub:
    """
    mkdir hisat2
    """
}
