process MAD_ROOTING {
    tag "${tree_file.baseName}"
    label 'process_low'
    
    publishDir "${params.outdir}/${unique_outdir}/rooted_trees", mode: params.publish_dir_mode

    input:
    tuple val(window_name), path(tree_file), path(alignment_file)
    val unique_outdir

    output:
    tuple val(window_name), path("${window_name}_rooted.tree"), path("${window_name}_rooting.log"), emit: rooted_trees
    path "*.log", emit: logs

    when:
    task.ext.when == null || task.ext.when

    script:    """
    # Set parameters with defaults
    ROOTING_STRATEGY="${params.rootdigger_strategy ?: 'modified-mad'}"
    EXHAUSTIVE_MODE="${params.rootdigger_exhaustive ?: false}"

    # Clean up taxa names to ensure consistency between alignment and tree
    # Extract only the accession part (before the pipe symbol) from alignment headers
    sed 's/^>\\([^|]*\\).*/\\>\\1/' ${alignment_file} > clean_alignment.fasta
    
    # Clean up tree file - replace full headers with just accession numbers
    sed 's/\\(LC[0-9]*\\)|[^:,()]*/\\1/g' ${tree_file} > clean_tree.newick
    
    # Run RootDigger for MAD-based rooting with configurable strategy
    if [ "\$EXHAUSTIVE_MODE" = "true" ]; then
        rootdigger \\
            --msa clean_alignment.fasta \\
            --tree clean_tree.newick \\
            --threads ${task.cpus} \\
            --initial-root-strategy "\$ROOTING_STRATEGY" \\
            --exhaustive \\
            --silent > rootdigger_output.txt
    else
        rootdigger \\
            --msa clean_alignment.fasta \\
            --tree clean_tree.newick \\
            --threads ${task.cpus} \\
            --initial-root-strategy "\$ROOTING_STRATEGY" \\
            --silent > rootdigger_output.txt
    fi

    # Extract just the tree (last line) from the output
    tail -1 rootdigger_output.txt > ${window_name}_rooted.tree

    # Check if rooting was successful
    if [ -s "${window_name}_rooted.tree" ]; then
        echo "Window: ${window_name}" > ${window_name}_rooting.log
        echo "RootDigger MAD rooting completed successfully" >> ${window_name}_rooting.log
        echo "Strategy: \$ROOTING_STRATEGY" >> ${window_name}_rooting.log
        echo "Exhaustive mode: \$EXHAUSTIVE_MODE" >> ${window_name}_rooting.log
        echo "Tree rooted successfully with RootDigger" >> ${window_name}_rooting.log
        
        # Log rooting confidence if available
        if grep -q "Log-Likelihood" rootdigger_output.txt; then
            echo "Root placement confidence information available" >> ${window_name}_rooting.log
            grep -E "(Log-Likelihood|Support)" rootdigger_output.txt >> ${window_name}_rooting.log
        fi
    else
        echo "RootDigger rooting failed for ${window_name}" >&2
        touch ${window_name}_rooted.tree
        echo "ERROR: RootDigger rooting failed" > ${window_name}_rooting.log
        exit 1
    fi
    """
}
