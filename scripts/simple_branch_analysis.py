#!/usr/bin/env python3
"""
Simple Branch Length Analysis and Correction Script

This script analyzes and fixes extremely long branch lengths in phylogenetic trees
that make visualization difficult.
"""

import re
import numpy as np
from pathlib import Path


def extract_branch_lengths(newick_string):
    """Extract all branch lengths from a Newick string."""
    pattern = r":([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
    matches = re.findall(pattern, newick_string)
    return [float(x) for x in matches]


def analyze_branch_lengths(tree_file):
    """Analyze branch lengths across all trees in a file."""
    print(f"Analyzing branch lengths in {tree_file}")

    all_lengths = []
    tree_count = 0
    problematic_trees = []

    with open(tree_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line:
                tree_count += 1
                lengths = extract_branch_lengths(line)
                all_lengths.extend(lengths)

                # Check for extremely long branches
                max_length = max(lengths) if lengths else 0
                if max_length > 1.0:  # Threshold for "problematic" branches
                    problematic_trees.append((line_num, max_length, lengths))

    # Statistical analysis
    all_lengths = np.array(all_lengths)

    print("\n=== BRANCH LENGTH ANALYSIS ===")
    print(f"Total trees analyzed: {tree_count}")
    print(f"Total branches analyzed: {len(all_lengths)}")
    print(f"Trees with branches > 1.0: {len(problematic_trees)}")

    print("\n=== BRANCH LENGTH STATISTICS ===")
    print(f"Mean: {np.mean(all_lengths):.6f}")
    print(f"Median: {np.median(all_lengths):.6f}")
    print(f"Min: {np.min(all_lengths):.6f}")
    print(f"Max: {np.max(all_lengths):.6f}")
    print(f"Standard deviation: {np.std(all_lengths):.6f}")

    # Percentiles
    percentiles = [90, 95, 99, 99.9]
    print("\n=== PERCENTILES ===")
    for p in percentiles:
        print(f"{p}th percentile: {np.percentile(all_lengths, p):.6f}")

    # Identify extreme outliers
    q99 = np.percentile(all_lengths, 99)
    extreme_outliers = all_lengths[all_lengths > q99 * 10]
    print(f"\nExtreme outliers (>10x 99th percentile): {len(extreme_outliers)}")
    if len(extreme_outliers) > 0:
        print(f"Extreme outlier values: {sorted(extreme_outliers, reverse=True)[:10]}")

    # Show problematic trees
    print("\n=== PROBLEMATIC TREES ===")
    for tree_num, max_len, lengths in problematic_trees[:5]:
        extreme_count = sum(1 for x in lengths if x > 1.0)
        print(
            f"Tree {tree_num}: Max branch = {max_len:.6f}, {extreme_count} branches > 1.0"
        )

    return all_lengths, problematic_trees


def identify_problematic_sequences(tree_file):
    """Identify which sequences are causing extremely long branches."""
    print("\n=== IDENTIFYING PROBLEMATIC SEQUENCES ===")

    sequence_max_branches = {}

    with open(tree_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line:
                # Extract sequence names and their associated branch lengths
                pattern = r"([A-Z0-9_\.]+):([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
                matches = re.findall(pattern, line)

                for seq_name, branch_len in matches:
                    branch_len = float(branch_len)
                    if seq_name not in sequence_max_branches:
                        sequence_max_branches[seq_name] = []
                    sequence_max_branches[seq_name].append(branch_len)

    # Calculate statistics for each sequence
    seq_stats = {}
    for seq_name, branches in sequence_max_branches.items():
        seq_stats[seq_name] = {
            "max": max(branches),
            "mean": np.mean(branches),
            "count": len(branches),
            "extreme_count": sum(1 for x in branches if x > 1.0),
        }

    # Sort by maximum branch length
    sorted_seqs = sorted(seq_stats.items(), key=lambda x: x[1]["max"], reverse=True)

    print("Top 10 sequences with longest branches:")
    for seq_name, stats in sorted_seqs[:10]:
        print(
            f"{seq_name}: max={stats['max']:.6f}, mean={stats['mean']:.6f}, "
            f"extreme_branches={stats['extreme_count']}/{stats['count']}"
        )

    return seq_stats


def cap_branch_lengths(tree_file, output_file, max_length=0.5):
    """Cap extremely long branch lengths to a maximum value."""
    print("\n=== CAPPING BRANCH LENGTHS ===")
    print(f"Maximum allowed branch length: {max_length}")

    capped_count = 0
    total_branches = 0

    with open(tree_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            line = line.strip()
            if line:
                # Find and replace long branch lengths
                def replace_branch(match):
                    nonlocal capped_count, total_branches
                    total_branches += 1
                    length = float(match.group(1))
                    if length > max_length:
                        capped_count += 1
                        return f":{max_length}"
                    return match.group(0)

                pattern = r":([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
                new_line = re.sub(pattern, replace_branch, line)
                outfile.write(new_line + "\n")

    print(f"Capped {capped_count} out of {total_branches} branches")
    print(f"Output written to: {output_file}")

    return capped_count, total_branches


def main():
    # Configuration
    tree_file = "results_200_15/best_rooted_trees.newick"

    if not Path(tree_file).exists():
        print(f"Error: Tree file {tree_file} not found!")
        return

    # Analyze branch lengths
    all_lengths, problematic_trees = analyze_branch_lengths(tree_file)

    # Identify problematic sequences
    seq_stats = identify_problematic_sequences(tree_file)

    # Create corrected versions
    print("\n=== CREATING CORRECTED TREE FILES ===")

    # Version 1: Cap at 0.5 (moderate correction)
    cap_branch_lengths(tree_file, "corrected_trees_capped_0.5.newick", max_length=0.5)

    # Version 2: Cap at 0.1 (aggressive correction for visualization)
    cap_branch_lengths(tree_file, "corrected_trees_capped_0.1.newick", max_length=0.1)

    # Version 3: Cap at 0.05 (very aggressive for clear visualization)
    cap_branch_lengths(tree_file, "corrected_trees_capped_0.05.newick", max_length=0.05)

    print("\n=== RECOMMENDATIONS ===")
    print("1. For basic visualization: Use corrected_trees_capped_0.5.newick")
    print("2. For clear visualization: Use corrected_trees_capped_0.1.newick")
    print("3. For very clear visualization: Use corrected_trees_capped_0.05.newick")
    print("4. The extreme branch lengths suggest:")
    print("   - Possible sequence quality issues")
    print("   - Alignment problems")
    print("   - Genuine biological divergence")
    print("   - Consider removing problematic sequences if appropriate")


if __name__ == "__main__":
    main()
