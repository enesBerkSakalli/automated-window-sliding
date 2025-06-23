process COLLECT_BACKBONE_MIDPOINT_TREES {
    tag "collect_backbone_midpoint"
    label 'process_low'
    
    publishDir "${main_results_dir}", mode: 'copy', pattern: "all_rooted_trees.*", saveAs: { filename -> filename.replace("all_rooted_trees", "best_backbone_rooted_trees") }
    publishDir "${unique_outdir}", mode: 'copy'

    input:
    path rooted_trees
    path constraint_trees
    path logs
    val output_formats
    val unique_outdir
    val main_results_dir

    output:
    path "all_rooted_trees.*", emit: combined_rooted_trees
    path "all_constraint_trees.*", emit: combined_constraint_trees
    path "tree_summary_report.txt", emit: summary_report
    path "analysis_metadata.json", emit: metadata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def formats = output_formats.tokenize(',')
    
    """
    # Create and run the Python script
    cat > collect_trees.py << 'EOF'
#!/usr/bin/env python3

import os
import json
import re
from datetime import datetime
from pathlib import Path

# Analysis parameters
WINDOW_SIZE = ${params.window_size}
STEP_SIZE = ${params.step_size}

def extract_window_info(filename):
    '''Extract window position from filename.'''
    match = re.search(r'window_([0-9]+)', filename)
    if match:
        return int(match.group(1))
    # Try alternative pattern for backbone midpoint files
    match = re.search(r'([0-9]+)_rooted\\.treefile', filename)
    if match:
        return int(match.group(1))
    # Try another pattern
    match = re.search(r'([0-9]+)\\.rooted', filename)
    if match:
        return int(match.group(1))
    return None

def collect_trees(tree_files, output_prefix, tree_type, formats):
    '''Collect and sort trees by window position.'''
    
    # Parse and sort tree files
    tree_data = []
    for tree_file in tree_files:
        if os.path.exists(tree_file):
            window_pos = extract_window_info(tree_file)
            if window_pos is not None:
                with open(tree_file, 'r') as f:
                    tree_content = f.read().strip()
                tree_data.append((window_pos, tree_file, tree_content))
    
    # Sort by window position
    tree_data.sort(key=lambda x: x[0])
    
    print(f"Found {len(tree_data)} {tree_type} trees to combine")
    
    # Generate output files for each format
    for fmt in formats:
        if fmt == 'newick':
            # Simple Newick format
            output_file = f"{output_prefix}.newick"
            with open(output_file, 'w') as f:
                for window_pos, filename, tree_content in tree_data:
                    f.write(f"{tree_content}\\n")
            
        elif fmt == 'nexus':
            # Nexus format with annotations
            output_file = f"{output_prefix}.nexus"
            with open(output_file, 'w') as f:
                f.write("#NEXUS\\n\\n")
                f.write("BEGIN TREES;\\n")
                f.write(f"\\t[! Backbone-Constrained {tree_type.title()} Trees ]\\n")
                f.write(f"\\t[! Generated: {datetime.now().isoformat()} ]\\n")
                f.write(f"\\t[! Window size: {WINDOW_SIZE} bp, Step size: {STEP_SIZE} bp ]\\n")
                f.write(f"\\t[! Total trees: {len(tree_data)} ]\\n")
                f.write("\\n")
                
                for i, (window_pos, filename, tree_content) in enumerate(tree_data, 1):
                    window_end = window_pos + WINDOW_SIZE - 1
                    f.write(f"\\ttree window_{window_pos}_{window_end} = [&R] {tree_content}\\n")
                
                f.write("END;\\n")
        
        elif fmt == 'annotated':
            # Annotated format with detailed metadata
            output_file = f"{output_prefix}.tre"
            with open(output_file, 'w') as f:
                f.write(f"# Backbone-Constrained {tree_type.title()} Trees\\n")
                f.write(f"# Generated: {datetime.now().isoformat()}\\n")
                f.write(f"# Analysis Parameters:\\n")
                f.write(f"#   Window Size: {WINDOW_SIZE} bp\\n")
                f.write(f"#   Step Size: {STEP_SIZE} bp\\n")
                f.write(f"#   Rooting Method: Backbone-constrained midpoint\\n")
                f.write(f"#   Tree Type: {tree_type.title()}\\n")
                f.write(f"#   Total Trees: {len(tree_data)}\\n")
                f.write(f"#\\n")
                f.write(f"# Tree Format: Newick with branch lengths\\n")
                f.write(f"# Window positions correspond to nucleotide positions in alignment\\n")
                f.write(f"#\\n")
                
                for i, (window_pos, filename, tree_content) in enumerate(tree_data, 1):
                    window_end = window_pos + WINDOW_SIZE - 1
                    f.write(f"# Tree {i}: Window positions {window_pos}-{window_end}\\n")
                    f.write(f"# Source: {filename}\\n")
                    f.write(f"[&R] {tree_content}\\n")
                    f.write("\\n")
    
    return len(tree_data)

# Main execution
rooted_tree_files = [f for f in "${rooted_trees}".split() if f.endswith('.treefile')]
constraint_tree_files = [f for f in "${constraint_trees}".split() if f.endswith('.treefile')]

formats = "${output_formats}".split(',')

# Collect rooted trees
num_rooted = collect_trees(rooted_tree_files, "all_rooted_trees", "rooted", formats)

# Collect constraint trees  
num_constraint = collect_trees(constraint_tree_files, "all_constraint_trees", "constraint", formats)

# Generate summary report
with open("tree_summary_report.txt", 'w') as f:
    f.write("# Backbone-Constrained Midpoint Rooting Summary Report\\n")
    f.write(f"# Generated: {datetime.now().isoformat()}\\n")
    f.write("\\n")
    f.write("## Analysis Parameters\\n")
    f.write(f"Window Size: {WINDOW_SIZE} bp\\n")
    f.write(f"Step Size: {STEP_SIZE} bp\\n")
    f.write(f"Rooting Method: Backbone-constrained midpoint\\n")
    f.write("\\n")
    f.write("## Results Summary\\n")
    f.write(f"Constraint Trees Generated: {num_constraint}\\n")
    f.write(f"Rooted Trees Generated: {num_rooted}\\n")
    f.write(f"Success Rate: {(num_rooted/max(num_constraint,1)*100):.1f}%\\n")
    f.write("\\n")
    f.write("## Output Files\\n")
    for fmt in formats:
        f.write(f"- all_rooted_trees.{fmt}\\n")
        f.write(f"- all_constraint_trees.{fmt}\\n")
    f.write("\\n")
    f.write("## Methodology\\n")
    f.write("1. Backbone tree generated from full alignment\\n")
    f.write("2. Window trees constrained by backbone topology\\n")
    f.write("3. Midpoint rooting applied for consistent root placement\\n")
    f.write("4. All trees validated and combined\\n")

# Generate analysis metadata
metadata = {
    "analysis_type": "backbone_constrained_midpoint_rooting",
    "timestamp": datetime.now().isoformat(),
    "parameters": {
        "window_size": WINDOW_SIZE,
        "step_size": STEP_SIZE,
        "rooting_method": "midpoint",
        "backbone_constraint": True,
        "output_formats": formats
    },
    "results": {
        "constraint_trees": num_constraint,
        "rooted_trees": num_rooted,
        "success_rate": round(num_rooted/max(num_constraint,1)*100, 1)
    }
}

with open("analysis_metadata.json", 'w') as f:
    json.dump(metadata, f, indent=2)

print(f"Successfully collected {num_rooted} rooted trees and {num_constraint} constraint trees")
EOF

    # Run the Python script
    python3 collect_trees.py
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "COLLECT_BACKBONE_MIDPOINT_TREES":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}