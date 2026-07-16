process CHROMAP_CHROMAP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed15f835e0883dfa27385f663655ad76280cba832c5a2db280ecd69f46877e91/data' :
        'community.wave.seqera.io/library/chromap_samtools:f82d647d01e30c96' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    path barcodes
    path whitelist
    path chr_order
    path pairs_chr_order

    output:
    tuple val(meta), path("*.bed.gz")     , optional:true, emit: bed
    tuple val(meta), path("*.bam")        , optional:true, emit: bam
    tuple val(meta), path("*.tagAlign.gz"), optional:true, emit: tagAlign
    tuple val(meta), path("*.pairs.gz")   , optional:true, emit: pairs
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    // Patch chromap/chromap module to include read group (otherwise picard/markduplicates erroring)
    // https://github.com/nf-core/chipseq/commit/b209c18d16e4b81dd0033014d8b76cbf077aa3d7#diff-bf809c928cb10b54d251dec8a140a2d5505150f97175d8d7b51dda9cc57971feR265
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()

    def file_extension = args.contains("--SAM") ? 'sam' : args.contains("--TagAlign")? 'tagAlign' : args.contains("--pairs")? 'pairs' : 'bed'
    if (barcodes) {
        args_list << "-b ${barcodes.join(',')}"
        if (whitelist) {
            args_list << "--barcode-whitelist $whitelist"
        }
    }
    if (chr_order) {
        args_list << "--chr-order $chr_order"
    }
    if (pairs_chr_order){
        args_list << "--pairs-natural-chr-order $pairs_chr_order"
    }
    def final_args = args_list.join(' ')
    def compression_cmds = "gzip -n ${prefix}.${file_extension}"
    if (args.contains("--SAM")) {
        compression_cmds = """
        samtools addreplacerg $args3 -o ${prefix}.rg.${file_extension} ${prefix}.${file_extension}
        samtools view $args2 -@ $task.cpus -bh \\
            -o ${prefix}.bam ${prefix}.rg.${file_extension}
        rm ${prefix}.rg.${file_extension}
        rm ${prefix}.${file_extension}
        """
    }
    if (meta.single_end) {
        """
        chromap \\
            $final_args \\
            -t $task.cpus \\
            -x $index \\
            -r $fasta \\
            -1 ${reads.join(',')} \\
            -o ${prefix}.${file_extension}
        
        $compression_cmds

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            chromap: \$(echo \$(chromap --version 2>&1))
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        chromap \\
            $final_args \\
            -t $task.cpus \\
            -x $index \\
            -r $fasta \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${prefix}.${file_extension}

        $compression_cmds

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            chromap: \$(echo \$(chromap --version 2>&1))
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.bed.gz
    touch ${prefix}.bam
    echo "" | gzip > ${prefix}.tagAlign.gz
    echo "" | gzip > ${prefix}.pairs.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromap: \$(echo \$(chromap --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
