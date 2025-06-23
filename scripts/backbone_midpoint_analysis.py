#!/usr/bin/env python3
"""
Optimized sliding window analysis with backbone constraint and midpoint rooting.
This version uses midpoint rooting instead of RootDigger for faster execution.
"""

import os
import subprocess
import sys
from Bio import Phylo
from io import StringIO


def run_iqtree_with_constraint(alignment_file, constraint_file, output_prefix):
    """
    Run IQ-TREE with backbone constraint and generate unrooted tree.

    Args:
        alignment_file: Path to alignment file
        constraint_file: Path to constraint tree file
        output_prefix: Prefix for output files

    Returns:
        bool: True if successful, False otherwise
    """
    cmd = [
        "iqtree2",
        "-s",
        alignment_file,
        "-g",
        constraint_file,
        "--prefix",
        output_prefix,
        "-nt",
        "AUTO",
        "-quiet",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        return result.returncode == 0
    except Exception as e:
        print(f"Error running IQ-TREE: {e}")
        return False


def midpoint_root_tree(tree_file, output_file):
    """
    Apply midpoint rooting to a tree using Biopython.

    Args:
        tree_file: Path to input unrooted tree
        output_file: Path to output rooted tree

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Read tree
        with open(tree_file, "r") as f:
            tree_str = f.read().strip()

        tree = Phylo.read(StringIO(tree_str), "newick")

        # Apply midpoint rooting
        tree.root_at_midpoint()

        # Write rooted tree
        with open(output_file, "w") as f:
            Phylo.write(tree, f, "newick")

        return True

    except Exception as e:
        print(f"Error in midpoint rooting: {e}")
        return False


def process_window_with_midpoint(
    alignment_file, constraint_file, window_id, output_dir
):
    """
    Process a single window with constraint-guided tree reconstruction and midpoint rooting.

    Args:
        alignment_file: Path to cleaned alignment file
        constraint_file: Path to backbone constraint tree
        window_id: Window identifier (e.g., "1", "21", etc.)
        output_dir: Directory for output files

    Returns:
        bool: True if successful, False otherwise
    """
    # Create output directories
    constraint_dir = os.path.join(output_dir, "constraint_trees_midpoint")
    rooted_dir = os.path.join(output_dir, "rooted_trees_midpoint")

    os.makedirs(constraint_dir, exist_ok=True)
    os.makedirs(rooted_dir, exist_ok=True)

    # Output files
    iqtree_prefix = os.path.join(constraint_dir, f"window_{window_id}_constrained")
    constraint_tree = f"{iqtree_prefix}.treefile"
    rooted_tree = os.path.join(
        rooted_dir, f"rooted_window_{window_id}_midpoint.treefile"
    )

    print(f"Processing window {window_id}...")

    # Step 1: Generate constraint-guided tree
    if not run_iqtree_with_constraint(alignment_file, constraint_file, iqtree_prefix):
        print(f"✗ IQ-TREE failed for window {window_id}")
        return False

    if not os.path.exists(constraint_tree):
        print(f"✗ Tree file not generated for window {window_id}")
        return False

    # Step 2: Apply midpoint rooting
    if not midpoint_root_tree(constraint_tree, rooted_tree):
        print(f"✗ Midpoint rooting failed for window {window_id}")
        return False

    print(f"✓ Window {window_id} completed successfully")
    return True


def main():
    """
    Run backbone-constrained analysis with midpoint rooting on all windows.
    """
    # Paths
    results_base = "results_optimized_200_20"
    new_results_dir = "results_backbone_midpoint_200_20"
    cleaned_alignments_dir = os.path.join(
        results_base, "sliding_windows", "alignments_cleaned"
    )
    backbone_tree = os.path.join(
        results_base, "backbone_analysis", "backbone_tree.treefile"
    )

    # Check prerequisites
    if not os.path.exists(cleaned_alignments_dir):
        print(
            f"Error: Cleaned alignments directory not found: {cleaned_alignments_dir}"
        )
        sys.exit(1)

    if not os.path.exists(backbone_tree):
        print(f"Error: Backbone tree not found: {backbone_tree}")
        sys.exit(1)

    # Create new results directory
    os.makedirs(new_results_dir, exist_ok=True)

    # Copy backbone tree for reference
    backbone_dest = os.path.join(new_results_dir, "backbone_tree.treefile")
    subprocess.run(["cp", backbone_tree, backbone_dest])

    # Find all cleaned alignment files
    alignment_files = []
    for file in os.listdir(cleaned_alignments_dir):
        if file.endswith(".fasta"):
            window_id = file.replace(".fasta", "")
            alignment_path = os.path.join(cleaned_alignments_dir, file)
            alignment_files.append((window_id, alignment_path))

    alignment_files.sort(key=lambda x: int(x[0]))
    print(f"Found {len(alignment_files)} alignment files to process")

    # Process each window
    successful_windows = 0
    failed_windows = 0

    for window_id, alignment_path in alignment_files:
        if process_window_with_midpoint(
            alignment_path, backbone_tree, window_id, new_results_dir
        ):
            successful_windows += 1
        else:
            failed_windows += 1

    # Summary
    print(f"\n{'=' * 60}")
    print("BACKBONE-CONSTRAINED MIDPOINT ROOTING SUMMARY")
    print(f"{'=' * 60}")
    print(f"Total windows processed: {len(alignment_files)}")
    print(f"Successfully processed: {successful_windows}")
    print(f"Failed: {failed_windows}")
    print(f"Success rate: {(successful_windows / len(alignment_files) * 100):.1f}%")
    print(f"\nResults saved to: {new_results_dir}")

    if successful_windows > 0:
        print("\nOutput directories:")
        print(f"  - Constraint trees: {new_results_dir}/constraint_trees_midpoint/")
        print(f"  - Midpoint rooted trees: {new_results_dir}/rooted_trees_midpoint/")


if __name__ == "__main__":
    main()
