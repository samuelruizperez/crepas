process REPLISEQ_RTNORMALIZE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6cebb5ebc440c01fc102ab494c0dccb249998735d2cca4ab517a8945a9428a65/data' :
        'community.wave.seqera.io/library/bioconductor-dnacopy_bioconductor-preprocesscore_r-argparse_r-base_pruned:635dd0ad085a07e6' }"
    // This container's R installation hardcodes /usr/bin/which somewhere in its base package
    // loading (independent of $PATH or exported shell functions -- verified with `singularity
    // --cleanenv`), which the minimal Wave-built image doesn't ship, so any library() call
    // that needs the `utils` base package fails outright. Bind the host's own /usr/bin/which
    // in at the same path as a workaround, when running under singularity/apptainer and the
    // host actually has one.
    containerOptions {
        def which_bin = file('/usr/bin/which')
        (workflow.containerEngine in ['singularity', 'apptainer'] && which_bin.exists())
            ? "--bind ${which_bin}:/usr/bin/which"
            : ''
    }

    input:
    tuple val(meta), path(counts), val(phases), val(breps)

    output:
    tuple val(meta), path("*.RT.raw.bedGraph")   , emit: raw
    tuple val(meta), path("*.RT.smooth.bedGraph"), emit: smooth
    tuple val(meta), path("*.qc.txt")            , emit: qc
    tuple val("${task.process}"), val('r-base'), eval("R --version | head -1 | sed 's/R version //; s/ .*//'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    repliseq_rtnormalize.R \\
        --counts ${counts} \\
        --phases ${phases.join(' ')} \\
        --breps ${breps.join(' ')} \\
        --outdir ./ \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.RT.raw.bedGraph
    touch ${prefix}.RT.smooth.bedGraph
    touch ${prefix}.qc.txt
    """
}
