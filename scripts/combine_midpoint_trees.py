#!/usr/bin/env python3
"""
Combine all backbone-constrained midpoint rooted trees into single files.
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


def combine_midpoint_trees():
    """
    Combine all midpoint rooted trees into a single file with metadata.
    """
    rooted_trees_dir = "results_backbone_midpoint_200_20/rooted_trees_midpoint"
    output_file = "ALL_BACKBONE_MIDPOINT_TREES.tre"

    if not os.path.exists(rooted_trees_dir):
        print(f"Error: Directory not found: {rooted_trees_dir}")
        return

    # Find all rooted tree files
    tree_files = []
    for file in os.listdir(rooted_trees_dir):
        if file.endswith(".treefile") and "midpoint" in file:
            tree_path = os.path.join(rooted_trees_dir, file)
            window_pos = extract_window_info(file)
            tree_files.append((window_pos, file, tree_path))

    # Sort by window position
    tree_files.sort(key=lambda x: x[0] if x[0] is not None else 0)

    print(f"Found {len(tree_files)} midpoint rooted trees to combine")

    # Create combined file with header
    with open(output_file, "w") as outfile:
        # Write header with metadata
        outfile.write(
            "# Backbone-Constrained Sliding Window Analysis - Midpoint Rooted Trees\n"
        )
        outfile.write(f"# Generated: {datetime.now().isoformat()}\n")
        outfile.write("# Analysis Parameters:\n")
        outfile.write("#   Window Size: 200 bp\n")
        outfile.write("#   Step Size: 20 bp\n")
        outfile.write("#   Rooting Method: Midpoint rooting\n")
        outfile.write("#   Constraint Optimization: Backbone-guided (IQ-TREE)\n")
        outfile.write(f"#   Total Windows: {len(tree_files)}\n")
        outfile.write("#   Taxa Count: 38 Norovirus sequences\n")
        outfile.write("#\n")
        outfile.write("# Tree Format: Newick with branch lengths\n")
        outfile.write("# Each tree is midpoint rooted and backbone-constrained\n")
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
            outfile.write("\n")

    print(f"All {len(tree_files)} midpoint rooted trees combined into: {output_file}")

    # Also create a simple Newick format file (trees only)
    simple_output = "ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre"
    with open(simple_output, "w") as outfile:
        for i, (window_pos, filename, tree_path) in enumerate(tree_files, 1):
            with open(tree_path, "r") as treefile:
                tree_content = treefile.read().strip()
            outfile.write(f"{tree_content}\n")

    print(f"Simple Newick format saved to: {simple_output}")

    # Create a summary table
    summary_file = "BACKBONE_MIDPOINT_TREE_SUMMARY.txt"
    with open(summary_file, "w") as outfile:
        outfile.write(
            "Window_Number\tStart_Position\tEnd_Position\tFilename\tRooting_Method\n"
        )
        for i, (window_pos, filename, tree_path) in enumerate(tree_files, 1):
            end_pos = window_pos + 199
            outfile.write(f"{i}\t{window_pos}\t{end_pos}\t{filename}\tMidpoint\n")

    print(f"Tree summary table saved to: {summary_file}")

    return output_file, simple_output, summary_file


def validate_midpoint_trees():
    """
    Quick validation of the midpoint rooted trees.
    """
    rooted_trees_dir = "results_backbone_midpoint_200_20/rooted_trees_midpoint"

    if not os.path.exists(rooted_trees_dir):
        print(f"Error: Directory not found: {rooted_trees_dir}")
        return

    tree_files = [f for f in os.listdir(rooted_trees_dir) if f.endswith(".treefile")]
    tree_files.sort()

    print(f"\nValidation of {len(tree_files)} midpoint rooted trees:")
    print("=" * 60)

    valid_trees = 0
    for tree_file in tree_files:
        tree_path = os.path.join(rooted_trees_dir, tree_file)
        try:
            with open(tree_path, "r") as f:
                tree_content = f.read().strip()

            # Basic validation - check if it's a valid Newick format
            if tree_content and tree_content.endswith(";") and "(" in tree_content:
                valid_trees += 1
                window_id = extract_window_info(tree_file)
                print(f"✓ Window {window_id}: Valid tree format")
            else:
                print(f"✗ {tree_file}: Invalid tree format")

        except Exception as e:
            print(f"✗ {tree_file}: Error reading file - {e}")

    print(f"\nValidation Summary:")
    print(
        f"Valid trees: {valid_trees}/{len(tree_files)} ({valid_trees / len(tree_files) * 100:.1f}%)"
    )


if __name__ == "__main__":
    print("Processing backbone-constrained midpoint rooted trees...")

    # Validate trees first
    validate_midpoint_trees()

    # Combine trees
    print("\nCombining trees into single files...")
    combine_midpoint_trees()

    print("\nBackbone-constrained midpoint analysis files ready!")
