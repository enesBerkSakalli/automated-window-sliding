process BACKBONE_MIDPOINT_ROOTING {
    tag "${window_name}"
    label 'process_medium'
    
    publishDir "${unique_outdir}/backbone_midpoint", mode: 'copy'

    input:
    tuple val(window_name), path(window_tree), path(window_alignment)
    path backbone_tree
    val unique_outdir

    output:
    tuple val(window_name), path("${window_name}_rooted.treefile"), emit: rooted_trees
    tuple val(window_name), path("${window_name}_constraint.treefile"), emit: constraint_trees
    path "${window_name}_rooting.log", emit: logs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Copy backbone tree as constraint
    cp ${backbone_tree} ${window_name}_constraint.treefile

    # Apply midpoint rooting
    backbone_midpoint_rooting.py ${window_tree} ${window_name}_rooted.treefile > ${window_name}_rooting.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "BACKBONE_MIDPOINT_ROOTING":
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${window_name}_rooted.treefile
    touch ${window_name}_constraint.treefile
    touch ${window_name}_rooting.log
    touch versions.yml
    """
}