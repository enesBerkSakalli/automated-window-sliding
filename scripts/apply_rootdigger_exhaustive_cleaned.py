#!/usr/bin/env python3
"""
Apply RootDigger exhaustive rooting to constraint trees using cleaned alignments.
This version uses the cleaned alignme      print("="*60)
    print("RootDigger Exhaustive Rooting Summary:")
    print("="*60)rint("="*60)
    print("RootDigger Exhaustive Rooting Summary:")
    print("="*60) to ensure taxa name consistency.
"""

import os
import subprocess
import sys
from pathlib import Path


def run_rootdigger_exhaustive(tree_file, alignment_file, output_dir):
    """
    Run RootDigger with exhaustive rooting on a single tree.

    Args:
        tree_file: Path to input tree file
        alignment_file: Path to corresponding alignment file (cleaned)
        output_dir: Directory for output files

    Returns:
        bool: True if successful, False otherwise
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Extract base name for output files
    tree_base = Path(tree_file).stem
    output_prefix = os.path.join(output_dir, f"rooted_{tree_base}")

    # RootDigger command with exhaustive rooting
    cmd = [
        "rootdigger",
        "--exhaustive",
        "--msa",
        alignment_file,
        "--tree",
        tree_file,
        "--prefix",
        output_prefix,
    ]

    print(f"Running RootDigger on {tree_base}...")
    print(f"Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 minute timeout per tree
        )

        if result.returncode == 0:
            print(f"✓ RootDigger completed successfully for {tree_base}")
            return True
        else:
            print(f"✗ RootDigger failed for {tree_base}")
            print(f"STDERR: {result.stderr}")
            return False

    except subprocess.TimeoutExpired:
        print(f"✗ RootDigger timed out for {tree_base}")
        return False
    except Exception as e:
        print(f"✗ Error running RootDigger for {tree_base}: {e}")
        return False


def get_corresponding_alignment(tree_file, alignments_dir):
    """
    Find the corresponding cleaned alignment file for a tree file.

    Args:
        tree_file: Path to tree file (e.g., window_21_constrained.treefile)
        alignments_dir: Directory containing cleaned alignment files

    Returns:
        str: Path to corresponding alignment file, or None if not found
    """
    # Extract window number from tree filename
    tree_name = Path(tree_file).stem
    if "window_" in tree_name:
        # Extract number from "window_21_constrained" -> "21"
        window_num = tree_name.split("_")[1]
        alignment_file = os.path.join(alignments_dir, f"{window_num}.fasta")

        if os.path.exists(alignment_file):
            return alignment_file
        else:
            print(
                f"Warning: Alignment file not found for {tree_name}: {alignment_file}"
            )
            return None
    else:
        print(f"Warning: Could not extract window number from {tree_name}")
        return None


def main():
    """
    Apply RootDigger exhaustive rooting to all constraint trees using cleaned alignments.
    """
    # Paths
    results_dir = "results_optimized_200_20"
    constraint_trees_dir = os.path.join(results_dir, "constraint_trees")
    cleaned_alignments_dir = os.path.join(
        results_dir, "sliding_windows", "alignments_cleaned"
    )
    rooted_trees_dir = os.path.join(results_dir, "rooted_trees_exhaustive")

    # Check if directories exist
    if not os.path.exists(constraint_trees_dir):
        print(f"Error: Constraint trees directory not found: {constraint_trees_dir}")
        sys.exit(1)

    if not os.path.exists(cleaned_alignments_dir):
        print(
            f"Error: Cleaned alignments directory not found: {cleaned_alignments_dir}"
        )
        sys.exit(1)

    # Create output directory
    os.makedirs(rooted_trees_dir, exist_ok=True)

    # Find all constraint tree files
    tree_files = []
    for file in os.listdir(constraint_trees_dir):
        if file.endswith(".treefile") and "constrained" in file:
            tree_files.append(os.path.join(constraint_trees_dir, file))

    tree_files.sort()
    print(f"Found {len(tree_files)} constraint trees to root")

    # Process each tree
    successful_roots = 0
    failed_roots = 0

    for tree_file in tree_files:
        tree_name = Path(tree_file).name
        print(f"\nProcessing {tree_name}...")

        # Find corresponding cleaned alignment
        alignment_file = get_corresponding_alignment(tree_file, cleaned_alignments_dir)

        if alignment_file is None:
            print(f"Skipping {tree_name} - no corresponding alignment found")
            failed_roots += 1
            continue

        # Run RootDigger
        if run_rootdigger_exhaustive(tree_file, alignment_file, rooted_trees_dir):
            successful_roots += 1
        else:
            failed_roots += 1

    # Summary
    print(f"\n{'=' * 60}")
    print("RootDigger Exhaustive Rooting Summary:")
    print("=" * 60)
    print(f"Total trees processed: {len(tree_files)}")
    print(f"Successfully rooted: {successful_roots}")
    print(f"Failed to root: {failed_roots}")
    print(f"Output directory: {rooted_trees_dir}")

    if successful_roots > 0:
        print(f"\nRooted trees can be found in: {rooted_trees_dir}")
        print("Each tree has multiple output files including:")
        print("  - .rooted_tree: The best rooted tree")
        print("  - .exhaustive: Log of exhaustive search")
        print("  - .log: Detailed RootDigger log")


if __name__ == "__main__":
    main()
