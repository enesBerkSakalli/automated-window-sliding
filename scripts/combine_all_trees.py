#!/usr/bin/env python3
"""
Combine all rooted trees into a single file with proper labels and metadata.
"""

import os
import re
from datetime import datetime


def extract_window_info(filename):
    """Extract window position from filename."""
    match = re.search(r"window_(\d+)", filename)
    if match:
        return int(match.group(1))
    return None


def combine_all_trees():
    """
    Combine all rooted trees into a single file with metadata.
    """
    rooted_trees_dir = "results_optimized_200_20/rooted_trees_exhaustive"
    output_file = "ALL_ROOTED_TREES_OPTIMIZED.tre"

    if not os.path.exists(rooted_trees_dir):
        print(f"Error: Directory not found: {rooted_trees_dir}")
        return

    # Find all rooted tree files
    tree_files = []
    for file in os.listdir(rooted_trees_dir):
        if file.endswith(".rooted.tree"):
            tree_path = os.path.join(rooted_trees_dir, file)
            window_pos = extract_window_info(file)
            tree_files.append((window_pos, file, tree_path))

    # Sort by window position
    tree_files.sort(key=lambda x: x[0] if x[0] is not None else 0)

    print(f"Found {len(tree_files)} rooted trees to combine")

    # Create combined file with header
    with open(output_file, "w") as outfile:
        # Write header with metadata
        outfile.write(
            "# Optimized Sliding Window Phylogenetic Analysis - All Rooted Trees\n"
        )
        outfile.write(f"# Generated: {datetime.now().isoformat()}\n")
        outfile.write("# Analysis Parameters:\n")
        outfile.write("#   Window Size: 200 bp\n")
        outfile.write("#   Step Size: 20 bp\n")
        outfile.write("#   Rooting Method: RootDigger exhaustive\n")
        outfile.write("#   Constraint Optimization: Backbone-guided\n")
        outfile.write(f"#   Total Windows: {len(tree_files)}\n")
        outfile.write("#   Taxa Count: 38 Norovirus sequences\n")
        outfile.write("#\n")
        outfile.write("# Tree Format: Newick with branch lengths\n")
        outfile.write("# Each tree is labeled with its window position and number\n")
        outfile.write("#\n")
        outfile.write(
            "# Window positions correspond to nucleotide positions in the refined alignment\n"
        )
        outfile.write(
            "# Window 1 = positions 1-200, Window 2 = positions 21-220, etc.\n"
        )
        outfile.write("#\n")

        # Write each tree with labels
        for i, (window_pos, filename, tree_path) in enumerate(tree_files, 1):
            window_num = i

            # Read tree content
            with open(tree_path, "r") as treefile:
                tree_content = treefile.read().strip()

            # Write tree with metadata
            outfile.write(
                f"# Tree {window_num}: Window positions {window_pos}-{window_pos + 199}\n"
            )
            outfile.write(f"# Source file: {filename}\n")
            outfile.write(f"[&R] {tree_content}\n")
            outfile.write(f"\n")

    print(f"All {len(tree_files)} rooted trees combined into: {output_file}")

    # Also create a simple Newick format file (trees only)
    simple_output = "ALL_ROOTED_TREES_SIMPLE.tre"
    with open(simple_output, "w") as outfile:
        for i, (window_pos, filename, tree_path) in enumerate(tree_files, 1):
            with open(tree_path, "r") as treefile:
                tree_content = treefile.read().strip()
            outfile.write(f"{tree_content}\n")

    print(f"Simple Newick format saved to: {simple_output}")

    # Create a summary table
    summary_file = "TREE_SUMMARY_TABLE.txt"
    with open(summary_file, "w") as outfile:
        outfile.write("Window_Number\tStart_Position\tEnd_Position\tFilename\n")
        for i, (window_pos, filename, tree_path) in enumerate(tree_files, 1):
            end_pos = window_pos + 199
            outfile.write(f"{i}\t{window_pos}\t{end_pos}\t{filename}\n")

    print(f"Tree summary table saved to: {summary_file}")

    return output_file, simple_output, summary_file


if __name__ == "__main__":
    combine_all_trees()
