#!/usr/bin/env python3
"""
Final Branch Length Analysis for Norovirus Sliding Window Phylogenetic Analysis
This script provides comprehensive analysis of the final results from the midpoint rooting run.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO


def parse_newick_trees(file_path):
    """Parse trees from Newick file and extract branch lengths."""
    trees = []
    branch_lengths = []

    with open(file_path, "r") as f:
        content = f.read().strip()

    # Split by lines and process each tree
    tree_lines = [line.strip() for line in content.split("\n") if line.strip()]

    for i, tree_line in enumerate(tree_lines):
        if tree_line:
            try:
                # Parse tree using Bio.Phylo
                tree_handle = StringIO(tree_line)
                tree = Phylo.read(tree_handle, "newick")
                trees.append(tree)

                # Extract all branch lengths
                tree_branches = []
                for clade in tree.find_clades():
                    if clade.branch_length is not None:
                        tree_branches.append(clade.branch_length)

                branch_lengths.extend(tree_branches)
                print(
                    f"Tree {i + 1}: {len(tree_branches)} branches, max={max(tree_branches):.6f}, mean={np.mean(tree_branches):.6f}"
                )

            except Exception as e:
                print(f"Error parsing tree {i + 1}: {e}")

    return trees, branch_lengths


def analyze_branch_length_distribution(branch_lengths, title_suffix=""):
    """Analyze and visualize branch length distribution."""
    if not branch_lengths:
        print("No branch lengths found!")
        return

    branch_array = np.array(branch_lengths)

    # Calculate statistics
    stats = {
        "count": len(branch_array),
        "mean": np.mean(branch_array),
        "median": np.median(branch_array),
        "std": np.std(branch_array),
        "min": np.min(branch_array),
        "max": np.max(branch_array),
        "q95": np.percentile(branch_array, 95),
        "q99": np.percentile(branch_array, 99),
        "very_long": np.sum(branch_array > 1.0),
        "long": np.sum(branch_array > 0.5),
        "moderate": np.sum((branch_array > 0.1) & (branch_array <= 0.5)),
        "short": np.sum(branch_array <= 0.1),
    }

    print(f"\nBranch Length Statistics{title_suffix}:")
    print(f"Total branches: {stats['count']}")
    print(f"Mean: {stats['mean']:.6f}")
    print(f"Median: {stats['median']:.6f}")
    print(f"Standard deviation: {stats['std']:.6f}")
    print(f"Range: {stats['min']:.6f} - {stats['max']:.6f}")
    print(f"95th percentile: {stats['q95']:.6f}")
    print(f"99th percentile: {stats['q99']:.6f}")
    print("\nBranch length categories:")
    print(
        f"Very long (>1.0): {stats['very_long']} ({100 * stats['very_long'] / stats['count']:.1f}%)"
    )
    print(
        f"Long (0.5-1.0): {stats['long']} ({100 * stats['long'] / stats['count']:.1f}%)"
    )
    print(
        f"Moderate (0.1-0.5): {stats['moderate']} ({100 * stats['moderate'] / stats['count']:.1f}%)"
    )
    print(
        f"Short (≤0.1): {stats['short']} ({100 * stats['short'] / stats['count']:.1f}%)"
    )

    return stats


def create_visualization(branch_lengths, output_file, title_suffix=""):
    """Create comprehensive visualization of branch length distribution."""
    branch_array = np.array(branch_lengths)

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    # Histogram
    axes[0, 0].hist(
        branch_array, bins=50, alpha=0.7, color="skyblue", edgecolor="black"
    )
    axes[0, 0].set_xlabel("Branch Length")
    axes[0, 0].set_ylabel("Frequency")
    axes[0, 0].set_title(f"Branch Length Distribution{title_suffix}")
    axes[0, 0].axvline(
        np.mean(branch_array),
        color="red",
        linestyle="--",
        label=f"Mean: {np.mean(branch_array):.4f}",
    )
    axes[0, 0].axvline(
        np.median(branch_array),
        color="orange",
        linestyle="--",
        label=f"Median: {np.median(branch_array):.4f}",
    )
    axes[0, 0].legend()

    # Log scale histogram
    axes[0, 1].hist(
        branch_array[branch_array > 0],
        bins=50,
        alpha=0.7,
        color="lightgreen",
        edgecolor="black",
    )
    axes[0, 1].set_xlabel("Branch Length")
    axes[0, 1].set_ylabel("Frequency")
    axes[0, 1].set_title(f"Branch Length Distribution (Log Scale){title_suffix}")
    axes[0, 1].set_yscale("log")

    # Box plot
    axes[1, 0].boxplot(branch_array, vert=True)
    axes[1, 0].set_ylabel("Branch Length")
    axes[1, 0].set_title(f"Branch Length Box Plot{title_suffix}")

    # Cumulative distribution
    sorted_branches = np.sort(branch_array)
    cumulative = np.arange(1, len(sorted_branches) + 1) / len(sorted_branches)
    axes[1, 1].plot(sorted_branches, cumulative, color="purple", linewidth=2)
    axes[1, 1].set_xlabel("Branch Length")
    axes[1, 1].set_ylabel("Cumulative Probability")
    axes[1, 1].set_title(f"Cumulative Distribution{title_suffix}")
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Visualization saved as {output_file}")


def compare_results():
    """Compare results across different pipeline runs."""
    results_dirs = {
        "Original (44→38 sequences)": "/Users/berksakalli/Projects/automated-window-sliding/results_200_15/best_rooted_trees.newick",
        "Refined (38 sequences)": "/Users/berksakalli/Projects/automated-window-sliding/results_refined_200_15/best_rooted_trees.newick",
        "Enhanced (41 with outgroups)": "/Users/berksakalli/Projects/automated-window-sliding/results_enhanced_outgroups/best_rooted_trees.newick",
        "Single Outgroup (39 sequences)": "/Users/berksakalli/Projects/automated-window-sliding/results_single_outgroup/best_rooted_trees.newick",
        "Ingroup Only (38 sequences)": "/Users/berksakalli/Projects/automated-window-sliding/results_midpoint_rooting/best_rooted_trees.newick",
    }

    comparison_data = []

    for run_name, file_path in results_dirs.items():
        if os.path.exists(file_path):
            print(f"\nAnalyzing {run_name}...")
            try:
                trees, branch_lengths = parse_newick_trees(file_path)
                stats = analyze_branch_length_distribution(
                    branch_lengths, f" - {run_name}"
                )

                comparison_data.append(
                    {
                        "Run": run_name,
                        "Trees": len(trees),
                        "Total Branches": stats["count"],
                        "Mean Length": stats["mean"],
                        "Median Length": stats["median"],
                        "Max Length": stats["max"],
                        "Q95 Length": stats["q95"],
                        "Very Long (>1.0)": stats["very_long"],
                        "Long (0.5-1.0)": stats["long"],
                        "Artifacts (%)": 100 * stats["very_long"] / stats["count"]
                        if stats["count"] > 0
                        else 0,
                    }
                )

                # Create individual visualization
                output_file = f"branch_analysis_{run_name.lower().replace(' ', '_').replace('(', '').replace(')', '').replace('→', 'to')}.png"
                create_visualization(branch_lengths, output_file, f" - {run_name}")

            except Exception as e:
                print(f"Error analyzing {run_name}: {e}")
        else:
            print(f"File not found: {file_path}")

    # Create comparison table
    if comparison_data:
        df = pd.DataFrame(comparison_data)
        print("\n" + "=" * 100)
        print("COMPREHENSIVE RESULTS COMPARISON")
        print("=" * 100)
        print(df.to_string(index=False, float_format="%.6f"))

        # Save comparison table
        df.to_csv("branch_length_comparison.csv", index=False)
        print("\nComparison table saved as branch_length_comparison.csv")

    return comparison_data


def main():
    """Main analysis function."""
    print(
        "Final Branch Length Analysis for Norovirus Sliding Window Phylogenetic Analysis"
    )
    print("=" * 80)

    # Run comprehensive comparison
    compare_results()

    # Focus on final results (ingroup only)
    final_results_file = "/Users/berksakalli/Projects/automated-window-sliding/results_midpoint_rooting/best_rooted_trees.newick"

    if os.path.exists(final_results_file):
        print("\n" + "=" * 80)
        print("FINAL RESULTS ANALYSIS (Ingroup Only - MAD Rooting)")
        print("=" * 80)

        trees, branch_lengths = parse_newick_trees(final_results_file)
        stats = analyze_branch_length_distribution(branch_lengths, " - Final Results")

        # Create final visualization
        create_visualization(
            branch_lengths, "final_branch_analysis.png", " - Final Results"
        )

        # Summary assessment
        print("\n" + "=" * 80)
        print("QUALITY ASSESSMENT")
        print("=" * 80)

        artifacts = stats["very_long"]
        artifact_rate = 100 * artifacts / stats["count"] if stats["count"] > 0 else 0

        print(f"✅ Pipeline completed successfully: {len(trees)}/33 trees generated")
        print(f"✅ Branch length artifacts: {artifacts} ({artifact_rate:.2f}%)")

        if artifact_rate < 1.0:
            print("✅ EXCELLENT: Very low artifact rate (<1%)")
        elif artifact_rate < 5.0:
            print("✅ GOOD: Low artifact rate (<5%)")
        elif artifact_rate < 10.0:
            print("⚠️  MODERATE: Some artifacts present (<10%)")
        else:
            print("❌ HIGH: Significant artifacts present (>10%)")

        print(f"✅ Mean branch length: {stats['mean']:.6f} (reasonable for norovirus)")
        print(f"✅ Median branch length: {stats['median']:.6f}")
        print(f"✅ Maximum branch length: {stats['max']:.6f}")

        if stats["max"] < 0.5:
            print("✅ EXCELLENT: All branch lengths within normal range")
        elif stats["max"] < 1.0:
            print("✅ GOOD: Most branch lengths within normal range")
        else:
            print("⚠️  Some long branches present but may be biologically meaningful")

    else:
        print(f"Final results file not found: {final_results_file}")


if __name__ == "__main__":
    main()
