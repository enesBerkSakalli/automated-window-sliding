process MAD_ROOTING {
    tag "${tree_file.baseName}"
    label 'process_low'
    
    publishDir "${params.outdir}/rooted_trees", mode: params.publish_dir_mode

    conda "bioconda::rootdigger"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rootdigger:1.0.3--h9ee0642_0' :
        'biocontainers/rootdigger:1.0.3--h9ee0642_0' }"

    input:
    tuple val(window_name), path(tree_file), path(alignment_file)
    val unique_outdir

    output:
    tuple val(window_name), path("${window_name}_rooted.tree"), emit: rooted_trees
    path("${window_name}_rooting.log"), emit: logs

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Debug: Print input files to verify they exist
    echo "=== MAD ROOTING DEBUG INFO ==="
    echo "Processing window: ${window_name}"
    echo "Tree file: ${tree_file}"
    echo "Alignment file: ${alignment_file}"
    echo "Current directory contents:"
    ls -la
    echo "Tree file exists: \$(test -f "${tree_file}" && echo "YES" || echo "NO")"
    echo "Alignment file exists: \$(test -f "${alignment_file}" && echo "YES" || echo "NO")"
    echo "=============================="

    # Verify input files exist
    if [ ! -f "${tree_file}" ]; then
        echo "ERROR: Tree file ${tree_file} not found" >&2
        exit 1
    fi
    
    if [ ! -f "${alignment_file}" ]; then
        echo "ERROR: Alignment file ${alignment_file} not found" >&2
        exit 1
    fi

    # Set parameters with defaults
    ROOTING_STRATEGY="${params.rootdigger_strategy ?: 'modified-mad'}"
    EXHAUSTIVE_MODE="${params.rootdigger_exhaustive ?: false}"

    # Clean up taxa names to ensure consistency between alignment and tree
    # Extract only the accession part (before the first space) from alignment headers
    sed 's/^>\\([^ ]*\\) .*/\\>\\1/' ${alignment_file} > clean_alignment.fasta
    
    # Tree file should already have clean accession numbers from IQ-TREE
    cp ${tree_file} clean_tree.newick
    
    # Run RootDigger for MAD-based rooting with configurable strategy
    if [ "\$EXHAUSTIVE_MODE" = "true" ]; then
        rootdigger \\
            --msa clean_alignment.fasta \\
            --tree clean_tree.newick \\
            --threads ${task.cpus} \\
            --initial-root-strategy "\$ROOTING_STRATEGY" \\
            --exhaustive \\
            --silent > rootdigger_output.txt 2>&1
    else
        rootdigger \\
            --msa clean_alignment.fasta \\
            --tree clean_tree.newick \\
            --threads ${task.cpus} \\
            --initial-root-strategy "\$ROOTING_STRATEGY" \\
            --silent > rootdigger_output.txt 2>&1
    fi

    # Check if rootdigger ran successfully
    if [ \$? -eq 0 ]; then
        # Extract just the tree (last line) from the output
        tail -1 rootdigger_output.txt > ${window_name}_rooted.tree
        
        echo "Window: ${window_name}" > ${window_name}_rooting.log
        echo "RootDigger MAD rooting completed successfully" >> ${window_name}_rooting.log
        echo "Strategy: \$ROOTING_STRATEGY" >> ${window_name}_rooting.log
        echo "Exhaustive mode: \$EXHAUSTIVE_MODE" >> ${window_name}_rooting.log
        
        # Log rooting confidence if available
        if grep -q "Log-Likelihood" rootdigger_output.txt; then
            echo "Root placement confidence information available" >> ${window_name}_rooting.log
            grep -E "(Log-Likelihood|Support)" rootdigger_output.txt >> ${window_name}_rooting.log
        fi
    else
        echo "RootDigger rooting failed for ${window_name}" >&2
        cat rootdigger_output.txt >&2
        
        # Create a dummy rooted tree (copy original)
        cp ${tree_file} ${window_name}_rooted.tree
        
        echo "ERROR: RootDigger rooting failed, using unrooted tree" > ${window_name}_rooting.log
        echo "Error output:" >> ${window_name}_rooting.log
        cat rootdigger_output.txt >> ${window_name}_rooting.log
    fi
    """
}
