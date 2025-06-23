#!/usr/bin/env python3
"""
Sliding Window Phylogenetic Analysis Results Viewer
Analyzes and visualizes results from the automated-window-sliding pipeline
"""

import os
import sys
import dendropy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def load_trees(tree_file):
    """Load phylogenetic trees from Newick format"""
    try:
        trees = dendropy.TreeList.get(path=tree_file, schema="newick")
        print(f"✅ Loaded {len(trees)} trees from {tree_file}")
        return trees
    except Exception as e:
        print(f"❌ Error loading trees: {e}")
        return None


def analyze_tree_statistics(trees, windows_log):
    """Analyze basic statistics for each tree"""
    stats = []

    # Load window information
    windows_df = pd.read_csv(windows_log, sep="\t")

    for i, tree in enumerate(trees):
        if i < len(windows_df):
            window_info = windows_df.iloc[i]

            # Calculate tree statistics
            tree_length = tree.length()
            max_depth = (
                tree.max_distance_from_root() if tree.max_distance_from_root() else 0
            )
            num_taxa = len(tree.leaf_nodes())

            stats.append(
                {
                    "window": i + 1,
                    "start_pos": window_info["start"],
                    "end_pos": window_info["end"],
                    "window_length": window_info["win_len"],
                    "tree_length": tree_length,
                    "max_depth": max_depth,
                    "num_taxa": num_taxa,
                    "avg_branch_length": tree_length / len(tree.nodes())
                    if len(tree.nodes()) > 0
                    else 0,
                }
            )

    return pd.DataFrame(stats)


def plot_tree_metrics(stats_df, output_dir):
    """Create plots of tree metrics across windows"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(
        "Sliding Window Phylogenetic Analysis Results", fontsize=16, fontweight="bold"
    )

    # Plot 1: Tree length across windows
    axes[0, 0].plot(
        stats_df["window"], stats_df["tree_length"], "b-o", linewidth=2, markersize=6
    )
    axes[0, 0].set_xlabel("Window Number")
    axes[0, 0].set_ylabel("Total Tree Length")
    axes[0, 0].set_title("Tree Length Across Sliding Windows")
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: Maximum depth (root-to-tip distance)
    axes[0, 1].plot(
        stats_df["window"], stats_df["max_depth"], "r-s", linewidth=2, markersize=6
    )
    axes[0, 1].set_xlabel("Window Number")
    axes[0, 1].set_ylabel("Maximum Root-to-Tip Distance")
    axes[0, 1].set_title("Evolutionary Depth Across Windows")
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: Average branch length
    axes[1, 0].plot(
        stats_df["window"],
        stats_df["avg_branch_length"],
        "g-^",
        linewidth=2,
        markersize=6,
    )
    axes[1, 0].set_xlabel("Window Number")
    axes[1, 0].set_ylabel("Average Branch Length")
    axes[1, 0].set_title("Average Branch Length Across Windows")
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: Genomic position coverage
    axes[1, 1].fill_between(
        range(1, len(stats_df) + 1),
        stats_df["start_pos"],
        stats_df["end_pos"],
        alpha=0.6,
        color="purple",
    )
    axes[1, 1].set_xlabel("Window Number")
    axes[1, 1].set_ylabel("Genomic Position")
    axes[1, 1].set_title("Genomic Coverage of Sliding Windows")
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    output_file = os.path.join(output_dir, "sliding_window_analysis.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✅ Analysis plot saved to: {output_file}")

    return output_file


def generate_summary_report(stats_df, output_dir):
    """Generate a summary statistics report"""
    report_file = os.path.join(output_dir, "tree_analysis_summary.txt")

    with open(report_file, "w") as f:
        f.write("SLIDING WINDOW PHYLOGENETIC ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"Total windows analyzed: {len(stats_df)}\n")
        f.write(
            f"Genomic region covered: {stats_df['start_pos'].min()} - {stats_df['end_pos'].max()} bp\n"
        )
        f.write(f"Average window size: {stats_df['window_length'].mean():.1f} bp\n\n")

        f.write("TREE STATISTICS:\n")
        f.write("-" * 20 + "\n")
        f.write(
            f"Average tree length: {stats_df['tree_length'].mean():.6f} ± {stats_df['tree_length'].std():.6f}\n"
        )
        f.write(
            f"Average max depth: {stats_df['max_depth'].mean():.6f} ± {stats_df['max_depth'].std():.6f}\n"
        )
        f.write(
            f"Average branch length: {stats_df['avg_branch_length'].mean():.6f} ± {stats_df['avg_branch_length'].std():.6f}\n"
        )
        f.write(f"Number of taxa per tree: {stats_df['num_taxa'].iloc[0]}\n\n")

        f.write("VARIATION ANALYSIS:\n")
        f.write("-" * 20 + "\n")
        f.write(
            f"Tree length CV: {(stats_df['tree_length'].std() / stats_df['tree_length'].mean() * 100):.2f}%\n"
        )
        f.write(
            f"Max depth CV: {(stats_df['max_depth'].std() / stats_df['max_depth'].mean() * 100):.2f}%\n"
        )
        f.write(
            f"Branch length CV: {(stats_df['avg_branch_length'].std() / stats_df['avg_branch_length'].mean() * 100):.2f}%\n\n"
        )

        # Identify windows with extreme values
        f.write("NOTABLE WINDOWS:\n")
        f.write("-" * 20 + "\n")
        max_tree_length_idx = stats_df["tree_length"].idxmax()
        min_tree_length_idx = stats_df["tree_length"].idxmin()
        max_depth_idx = stats_df["max_depth"].idxmax()

        f.write(
            f"Longest tree: Window {stats_df.loc[max_tree_length_idx, 'window']} "
            f"(length: {stats_df.loc[max_tree_length_idx, 'tree_length']:.6f})\n"
        )
        f.write(
            f"Shortest tree: Window {stats_df.loc[min_tree_length_idx, 'window']} "
            f"(length: {stats_df.loc[min_tree_length_idx, 'tree_length']:.6f})\n"
        )
        f.write(
            f"Deepest tree: Window {stats_df.loc[max_depth_idx, 'window']} "
            f"(depth: {stats_df.loc[max_depth_idx, 'max_depth']:.6f})\n"
        )

    print(f"✅ Summary report saved to: {report_file}")
    return report_file


def main():
    """Main analysis function"""
    # Find the most recent results directory
    results_base = "/Users/berksakalli/Projects/automated-window-sliding/results"

    if not os.path.exists(results_base):
        print(f"❌ Results directory not found: {results_base}")
        return

    # Find the most recent analysis
    result_dirs = [
        d for d in os.listdir(results_base) if d.startswith("norovirus_sliding_window")
    ]
    if not result_dirs:
        print("❌ No analysis results found")
        return

    latest_result = sorted(result_dirs)[-1]
    result_path = os.path.join(results_base, latest_result)

    print(f"🔍 Analyzing results from: {latest_result}")

    # File paths
    rooted_trees_file = os.path.join(result_path, "best_rooted_trees.newick")
    unrooted_trees_file = os.path.join(result_path, "best_trees.newick")

    # Find windows.log file
    windows_log = None
    for root, dirs, files in os.walk(result_path):
        if "windows.log" in files:
            windows_log = os.path.join(root, "windows.log")
            break

    if not windows_log:
        print("❌ windows.log file not found")
        return

    # Load and analyze trees
    print("\n📊 Loading phylogenetic trees...")
    trees = load_trees(rooted_trees_file)

    if trees is None:
        return

    # Analyze tree statistics
    print("📈 Calculating tree statistics...")
    stats_df = analyze_tree_statistics(trees, windows_log)

    # Create output directory for analysis
    analysis_output_dir = "/Users/berksakalli/Projects/automated-window-sliding"

    # Generate plots
    print("📊 Creating analysis plots...")
    plot_file = plot_tree_metrics(stats_df, analysis_output_dir)

    # Generate summary report
    print("📝 Generating summary report...")
    report_file = generate_summary_report(stats_df, analysis_output_dir)

    # Save statistics as CSV
    csv_file = os.path.join(analysis_output_dir, "sliding_window_tree_stats.csv")
    stats_df.to_csv(csv_file, index=False)
    print(f"✅ Tree statistics saved to: {csv_file}")

    print(f"\n🎉 Analysis complete! Generated files:")
    print(f"   📊 Plot: {plot_file}")
    print(f"   📝 Report: {report_file}")
    print(f"   📋 Data: {csv_file}")

    # Display quick summary
    print(f"\n📋 Quick Summary:")
    print(f"   • {len(trees)} sliding windows analyzed")
    print(f"   • Average tree length: {stats_df['tree_length'].mean():.6f}")
    print(
        f"   • Tree length range: {stats_df['tree_length'].min():.6f} - {stats_df['tree_length'].max():.6f}"
    )
    print(
        f"   • Genomic coverage: {stats_df['start_pos'].min()}-{stats_df['end_pos'].max()} bp"
    )


if __name__ == "__main__":
    main()
