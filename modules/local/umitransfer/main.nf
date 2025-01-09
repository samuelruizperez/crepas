process UMITRANSFER {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi-transfer:1.5.0--h715e4b3_0' :
        'quay.io/biocontainers/umi-transfer:1.5.0--h715e4b3_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        umi-transfer \\
            external \\
            --in ${reads[0]} \\
            --in2 ${reads[0]} \\
            --umi ${reads[1]} \\
            --threads $task.cpus \\
            --out ${prefix}.umi_extract.fastq.gz \\
            --out2 '/dev/null' \\
            --gzip \\
            $args \\
            > ${prefix}.umi_extract.log
            
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umi-transfer external: \$( umi-transfer external --version | sed '/version:/!d; s/.*: //' )
        END_VERSIONS
        """
    }  else {
        """
        umi-transfer \\
            external \\
            --in ${reads[0]} \\
            --in2 ${reads[1]} \\
            --umi ${reads[2]} \\
            --threads $task.cpus \\
            --out ${prefix}.umi_extract_1.fastq.gz \\
            --out2 ${prefix}.umi_extract_2.fastq.gz \\
            --gzip \\
            $args \\
            > ${prefix}.umi_extract.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umi-transfer external: \$( umi-transfer external --version | sed '/version:/!d; s/.*: //' )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}.umi_extract.fastq.gz"
    } else {
        output_command = "echo '' | gzip > ${prefix}.umi_extract_1.fastq.gz ;"
        output_command += "echo '' | gzip > ${prefix}.umi_extract_2.fastq.gz"
    }
    """
    touch ${prefix}.umi_extract.log
    ${output_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi-transfer external: \$( umi-transfer external --version | sed '/version:/!d; s/.*: //' )
    END_VERSIONS
    """
}
