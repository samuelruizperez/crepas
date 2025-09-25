//
// Generating self-pseudoreplicates from TAGALIGN files
// Based on: https://github.com/ENCODE-DCC/chip-seq-pipeline2; step 2c
//
process TAGALIGN_SELF_PSEUDOREPLICATES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(tagalign)

    output:
    tuple val(meta), path("*spr1.tagAlign"), emit: tagalign1
    tuple val(meta), path("*spr2.tagAlign"), emit: tagalign2
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
        # Get total number of read pairs
        nlines=$( wc -l < ${tagalign} )
        nlines=$(( (nlines + 1) / 2 ))

        # Shuffle and split tagAlign file into 2 equal parts
        # Will produce prefix.shuf.split.tagAlign00 and prefix.shuf.split.tagAlign01
        shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(wc -c < ${tagalign}) -nosalt </dev/zero 2>/dev/null) ${tagalign} \\
            | split -d -l \${nlines} - ${prefix}.shuf.split.tagAlign

        # Convert reads into standard tagAlign file
        gzip -nc “${prefix}.shuf.split.tagAlign00" > ${prefix}.spr1.tagAlign
        rm "${prefix}.shuf.split.tagAlign00"
        gzip -nc “${prefix}.shuf.split.tagAlign01" > ${prefix}.spr2.tagAlign
        rm "${prefix}.shuf.split.tagAlign01"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
        END_VERSIONS
        """
    } else {
        """
        # Make temporary fake BEDPE file from FINAL_TA_FILE
        sed 'N;s/\\n/\\t/' ${tagalign} > ${prefix}.tmp.bedpe

        # Get total number of read pairs
        nlines=$( wc -l < ${prefix}.tmp.bedpe )
        nlines=$(( (nlines + 1) / 2 ))

        # Shuffle and split BEDPE file into 2 equal parts
        shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(wc -c < ${tagalign}) -nosalt </dev/zero 2>/dev/null) ${prefix}.tmp.bedpe \\
            | split -d -l \${nlines} - ${prefix}.shuf.split.tagAlign

        # Convert reads into standard tagAlign file
        awk 'BEGIN{OFS="\\t"}{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12} "${prefix}.shuf.split.tagAlign00" \\
            | gzip -nc > ${prefix}.spr1.tagAlign

        awk 'BEGIN{OFS="\\t"}{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12} "${prefix}.shuf.split.tagAlign01" \\
            | gzip -nc > ${prefix}.spr2.tagAlign
        
        rm "${prefix}.shuf.split.tagAlign00"
        rm "${prefix}.shuf.split.tagAlign01"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.spr1.tagAlign
    touch  ${prefix}.spr2.tagAlign

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
