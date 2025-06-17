process COLLECT_ROOTED_TREES {
    label 'process_single'
    
    publishDir "${params.outdir}/${unique_outdir}", mode: params.publish_dir_mode
    publishDir "${params.outdir}/", mode: params.publish_dir_mode, pattern: "rooted_trees_collection.*", saveAs: { filename -> filename.replace("rooted_trees_collection", "best_rooted_trees") }

    input:
    path rooted_trees
    path mad_logs
    val unique_outdir
    val output_format

    output:
    path "rooted_trees_collection.*", emit: rooted_collection
    path "mad_rooting_summary.txt", emit: rooting_summary
    path "root_quality_report.txt", emit: quality_report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    
    import os
    import re
    from pathlib import Path
    
    # Collect all rooted trees
    rooted_tree_files = []
    for tree_file in Path('.').glob('*_rooted.tree'):
        if tree_file.stat().st_size > 0:  # Only non-empty files
            rooted_tree_files.append(tree_file)
    
    rooted_tree_files.sort(key=lambda x: int(re.search(r'(\\d+)', x.stem).group(1)))
    
    # Output formats
    output_formats = "${output_format}".split(',')
    
    # Create collections based on requested formats
    if 'newick' in output_formats:
        with open('rooted_trees_collection.newick', 'w') as f:
            for tree_file in rooted_tree_files:
                with open(tree_file, 'r') as tf:
                    tree_content = tf.read().strip()
                    if tree_content and not tree_content.startswith('ERROR'):
                        f.write(tree_content + '\\n')
    
    if 'nexus' in output_formats:
        with open('rooted_trees_collection.nexus', 'w') as f:
            f.write('#NEXUS\\n\\n')
            f.write('BEGIN TREES;\\n')
            for i, tree_file in enumerate(rooted_tree_files, 1):
                with open(tree_file, 'r') as tf:
                    tree_content = tf.read().strip()
                    if tree_content and not tree_content.startswith('ERROR'):
                        window_name = tree_file.stem.replace('_rooted', '')
                        f.write(f'    TREE tree_{i} = [&R] {tree_content}\\n')
            f.write('END;\\n')
    
    # Generate RootDigger rooting summary
    with open('mad_rooting_summary.txt', 'w') as f:
        f.write('RootDigger MAD Rooting Summary\\n')
        f.write('=============================\\n\\n')
        f.write(f'Total trees processed: {len(rooted_tree_files)}\\n')
        f.write(f'Successfully rooted: {len([t for t in rooted_tree_files if Path(t).stat().st_size > 0])}\\n')
        f.write('\\nWindow-by-window results:\\n')
        f.write('Window\\tStatus\\tRooting_Strategy\\n')
        
        successful_roots = 0
        for tree_file in rooted_tree_files:
            window_name = tree_file.stem.replace('_rooted', '')
            log_file = tree_file.parent / f'{window_name}_rooting.log'
            
            if tree_file.stat().st_size > 0:
                status = 'SUCCESS'
                successful_roots += 1
                
                # Try to extract rooting strategy
                rooting_strategy = 'modified-mad'
                if log_file.exists():
                    with open(log_file, 'r') as lf:
                        log_content = lf.read()
                        # Look for strategy information
                        if 'Strategy:' in log_content:
                            for line in log_content.split('\\n'):
                                if 'Strategy:' in line:
                                    rooting_strategy = line.split('Strategy:')[1].strip()
                                    break
            else:
                status = 'FAILED'
                rooting_strategy = 'N/A'
            
            f.write(f'{window_name}\\t{status}\\t{rooting_strategy}\\n')
    
    # Generate quality report
    with open('root_quality_report.txt', 'w') as f:
        f.write('Root Quality Assessment Report\\n')
        f.write('=============================\\n\\n')
        f.write('This report provides quality metrics for MAD rooting.\\n')
        f.write('MAD (Minimal Ancestor Deviation) rooting is particularly suitable for:\\n')
        f.write('- Sliding window analysis without outgroups\\n')
        f.write('- Maintaining consistency across windows\\n')
        f.write('- Datasets with rate heterogeneity\\n\\n')
        
        f.write('Rooting Method: MAD (Minimal Ancestor Deviation)\\n')
        f.write(f'Windows successfully rooted: {successful_roots}/{len(rooted_tree_files)}\\n')
        f.write(f'Success rate: {(successful_roots/len(rooted_tree_files)*100):.1f}%\\n\\n')
        
        f.write('Quality Guidelines:\\n')
        f.write('- MAD rooting is deterministic and reproducible\\n')
        f.write('- Robust to rate heterogeneity across lineages\\n')
        f.write('- Suitable for phylogenomic sliding window analysis\\n')
        f.write('- For validation, consider RootDigger for uncertain cases\\n\\n')
        
        if successful_roots < len(rooted_tree_files):
            f.write('Failed windows may indicate:\\n')
            f.write('- Very short branches or star-like topologies\\n')
            f.write('- Insufficient phylogenetic signal\\n')
            f.write('- Consider concatenating failed windows or using alternative rooting\\n')
    
    print(f"MAD rooting completed: {successful_roots}/{len(rooted_tree_files)} trees successfully rooted")
    """
}
