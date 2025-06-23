process GENERATE_BACKBONE_TREE {
    tag "backbone_tree"
    label 'process_medium'
    
    publishDir "${unique_outdir}/backbone_analysis", mode: 'copy'

    input:
    path alignment
    val unique_outdir

    output:
    path "backbone_tree.treefile", emit: backbone_tree
    path "backbone_tree.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_backbone_tree.py "${alignment}" "${params.model ?: 'MFP'}" "${task.cpus}"
    """

    stub:
    """
    touch backbone_tree.treefile
    touch backbone_tree.log
    touch versions.yml
    """
}
