process DANPOS2_DPEAK {
    tag "${meta.id}"
    label 'process_low_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/dd9856b8473bf4b2e67285cc76dda3f9184e9bfe8432c8b667cde9ce1f971ef3/data'
        : 'community.wave.seqera.io/library/danpos:2.2.2--34b9ae5a26a3f9b4'}"

    input:
    tuple val(meta), 
          path(treatment,       stageAs: "treatment/*"), 
          path(treatment_input, stageAs: "treatment_input/*"), 
          val(treatment_count), 
          path(control,         stageAs: 'control/*'), 
          path(control_input,   stageAs: 'control_input/*'), 
          val(control_count)

    output:
    tuple val(meta), path("result/pooled/*.peaks.xls"), emit: pooled_xls
    tuple val(meta), path("result/*.peaks.integrative.xls"), emit: integrative_peaks, optional: true
    tuple val(meta), path("result/pooled/*input.wig"), emit: pooled_input_wig
    tuple val(meta), path("result/pooled/*smooth.wig"), emit: pooled_treat_wig
    tuple val(meta), path("result/pooled/*refregions.xls"), emit: pooled_bed, optional: true
    tuple val(meta), path("result/diff/*.wig"), emit: diff_wig, optional: true
    tuple val(meta), path("result/diff/*local_gain.refpeaks.xls"), emit: local_gain, optional: true
    tuple val(meta), path("result/diff/*local_loss.refpeaks.xls"), emit: local_loss, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_arg = control ? "treatment/:control/" : "treatment/"
    def treatment_input_arg = treatment_input ? "treatment/:treatment_input/" : ""
    def control_input_arg = control_input ? "control/:control_input/" : ""
    def input_arg = treatment_input || control_input ? "--bg ${treatment_input_arg + control_input_arg}" : ""
    def treatment_count_arg = treatment_count ? "treatment/${treatment_count}" : ""
    def control_count_arg = control_count ? "control/${control_count}" : ""
    def count_arg = treatment_count || control_count ? "--count ${treatment_count_arg + control_count_arg}" : ""
    def se = meta.single_end ? "" : "--paired 1"
    """
    danpos.py dpeak \\
        ${file_arg} \\
        ${args} \\
        ${input_arg} \\
        ${count_arg} \\
        ${se} \\
        --out result/

    # replace any occurrence of 'treatment' in filenames with the prefix
    find ./ -type f -name '*treatment*' -exec bash -c 'f="{}"; mv "\$f" "\${f//treatment/${prefix}}"' \\;

    # replace any occurrence of 'treatment' in subdirectory names with the prefix
    #find ./ -mindepth 1 -maxdepth 1 -type d -name '*treatment*' -exec bash -c 'd="{}"; mv "\$d" "\${d//treatment/${prefix}}"' \\;
    
    # TODO: this is just to keep track of the final dir structure:
    find .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DANPOS: \$(echo \$(danpos.py dpeak -h 2>&1) | grep -oP 'danpos\\s+\\K[\\d.]+' | head -1)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.peaks.xls
    touch ${prefix}.input.wig
    touch ${prefix}.smooth.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DANPOS: \$(echo \$(danpos.py dpeak -h 2>&1) | grep -oP 'danpos\\s+\\K[\\d.]+' | head -1)
    END_VERSIONS
    """
}
