process IGV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3':
        'biocontainers/python:3.8.3' }"

    input:
    tuple path(files), val(outpaths), val(colors), stageAs: "?/*"
    tuple path(fasta), val(fasta_outpath)

    output:
    path "*files.txt", emit: txt
    path "*.xml", emit: xml
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "igv_session."
    def outpath_list = outpaths.join(";")
    def color_list = colors.join(";")
    """
    # save outpaths and colors to a tab-separated file with one line per file
    paste \
        <(echo "${outpath_list}" | tr ';' '\\n') \
        <(echo "${color_list}" | tr ';' '\\n') \
        > ${prefix}.files.txt

    igv_files_to_session.py \
        ${args} \
        --file_list ${prefix}.files.txt \
        --genome_fasta ${fasta_outpath} \
        --xml_output ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "igv_session."
    """
    touch ${prefix}.files.txt
    touch ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
