#!/usr/bin/env python3
"""
Compare Topological Congruence: Original vs Improved
Analyzes and compares the Robinson-Foulds distances between approaches.
"""

import os
import numpy as np
import dendropy
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt
import pandas as pd


def load_and_analyze_trees(file_path, approach_name):
    """Load trees and calculate RF distances."""

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None, None

    print(f"\nAnalyzing {approach_name}...")
    print(f"Loading trees from {file_path}")

    trees = []
    with open(file_path, "r") as f:
        content = f.read().strip()

    tree_lines = [line.strip() for line in content.split("\n") if line.strip()]

    for i, tree_line in enumerate(tree_lines):
        if tree_line:
            try:
                tree = dendropy.Tree.get(data=tree_line, schema="newick")
                trees.append(tree)
            except Exception as e:
                print(f"Error loading tree {i + 1}: {e}")

    print(f"Loaded {len(trees)} trees")

    # Calculate RF distances
    n_trees = len(trees)
    adjacent_distances = []

    if n_trees < 2:
        return trees, []

    # Ensure all trees have the same taxon namespace
    taxon_namespace = trees[0].taxon_namespace
    for tree in trees[1:]:
        tree.migrate_taxon_namespace(taxon_namespace)

    print("Calculating adjacent RF distances...")
    for i in range(n_trees - 1):
        try:
            rf_distance = treecompare.robinson_foulds_distance(trees[i], trees[i + 1])
            adjacent_distances.append(rf_distance)
            if i < 10:  # Show first 10
                print(f"  Window {i + 1} → {i + 2}: RF = {rf_distance:.3f}")
            elif i == 10:
                print("  ...")
        except Exception as e:
            print(f"Error calculating RF between trees {i + 1} and {i + 2}: {e}")

    if adjacent_distances:
        print(f"\nAdjacent RF distances statistics:")
        print(f"  Count: {len(adjacent_distances)}")
        print(f"  Mean: {np.mean(adjacent_distances):.3f}")
        print(f"  Median: {np.median(adjacent_distances):.3f}")
        print(f"  Min: {np.min(adjacent_distances):.3f}")
        print(f"  Max: {np.max(adjacent_distances):.3f}")
        print(f"  Std: {np.std(adjacent_distances):.3f}")

    return trees, adjacent_distances


def compare_approaches():
    """Compare original and improved approaches."""

    print("Topological Congruence Comparison")
    print("=" * 50)

    approaches = {
        "Original (200bp/15bp)": {
            "file": "results_midpoint_rooting/best_rooted_trees.newick",
            "window_size": 200,
            "step_size": 15,
        },
        "Improved (400bp/5bp)": {
            "file": "results_improved_strategy2/best_rooted_trees.newick",
            "window_size": 400,
            "step_size": 5,
        },
    }

    results = {}

    for approach_name, params in approaches.items():
        trees, distances = load_and_analyze_trees(params["file"], approach_name)

        if distances:
            overlap_pct = (
                (params["window_size"] - params["step_size"]) / params["window_size"]
            ) * 100

            results[approach_name] = {
                "n_trees": len(trees),
                "n_distances": len(distances),
                "mean_rf": np.mean(distances),
                "median_rf": np.median(distances),
                "min_rf": np.min(distances),
                "max_rf": np.max(distances),
                "std_rf": np.std(distances),
                "window_size": params["window_size"],
                "step_size": params["step_size"],
                "overlap_pct": overlap_pct,
                "distances": distances,
            }

    return results


def create_comparison_visualization(
    results, output_file="topological_congruence_comparison.png"
):
    """Create comparison visualization."""

    if len(results) < 2:
        print("Need at least 2 approaches for comparison")
        return

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    colors = ["blue", "red", "green", "purple"]
    approach_names = list(results.keys())

    # 1. RF distance distributions
    for i, (approach, data) in enumerate(results.items()):
        if data["distances"]:
            axes[0, 0].hist(
                data["distances"],
                bins=20,
                alpha=0.7,
                label=f"{approach}\n(mean: {data['mean_rf']:.2f})",
                color=colors[i],
            )

    axes[0, 0].set_title("Robinson-Foulds Distance Distributions")
    axes[0, 0].set_xlabel("RF Distance")
    axes[0, 0].set_ylabel("Frequency")
    axes[0, 0].legend()

    # 2. Box plot comparison
    box_data = [data["distances"] for data in results.values()]
    box_labels = [name.split("(")[0].strip() for name in results.keys()]

    axes[0, 1].boxplot(box_data, labels=box_labels)
    axes[0, 1].set_title("RF Distance Box Plot Comparison")
    axes[0, 1].set_ylabel("RF Distance")
    axes[0, 1].tick_params(axis="x", rotation=45)

    # 3. Adjacent window RF distances over windows
    for i, (approach, data) in enumerate(results.items()):
        if data["distances"]:
            x_values = range(1, len(data["distances"]) + 1)
            axes[1, 0].plot(
                x_values,
                data["distances"],
                "o-",
                label=approach,
                color=colors[i],
                alpha=0.7,
            )

    axes[1, 0].set_title("Adjacent Window RF Distances")
    axes[1, 0].set_xlabel("Window Transition")
    axes[1, 0].set_ylabel("RF Distance")
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # 4. Summary statistics comparison
    metrics = ["mean_rf", "median_rf", "max_rf", "std_rf"]
    metric_labels = ["Mean RF", "Median RF", "Max RF", "Std RF"]

    x_pos = np.arange(len(metrics))
    width = 0.35

    for i, (approach, data) in enumerate(results.items()):
        values = [data[metric] for metric in metrics]
        axes[1, 1].bar(
            x_pos + i * width, values, width, label=approach, color=colors[i], alpha=0.7
        )

    axes[1, 1].set_title("Summary Statistics Comparison")
    axes[1, 1].set_xlabel("Metrics")
    axes[1, 1].set_ylabel("RF Distance")
    axes[1, 1].set_xticks(x_pos + width / 2)
    axes[1, 1].set_xticklabels(metric_labels)
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Comparison visualization saved as {output_file}")


def create_comparison_table(results):
    """Create comparison table."""

    table_data = []

    for approach, data in results.items():
        table_data.append(
            {
                "Approach": approach,
                "Trees": data["n_trees"],
                "Window Size": f"{data['window_size']}bp",
                "Step Size": f"{data['step_size']}bp",
                "Overlap": f"{data['overlap_pct']:.1f}%",
                "Mean RF": f"{data['mean_rf']:.3f}",
                "Median RF": f"{data['median_rf']:.3f}",
                "Min RF": f"{data['min_rf']:.3f}",
                "Max RF": f"{data['max_rf']:.3f}",
                "Std RF": f"{data['std_rf']:.3f}",
            }
        )

    df = pd.DataFrame(table_data)

    print("\n" + "=" * 120)
    print("TOPOLOGICAL CONGRUENCE COMPARISON TABLE")
    print("=" * 120)
    print(df.to_string(index=False))

    # Save to CSV
    df.to_csv("topological_congruence_comparison.csv", index=False)
    print(f"\nComparison table saved to topological_congruence_comparison.csv")

    return df


def assess_improvement(results):
    """Assess the improvement between approaches."""

    if len(results) < 2:
        return

    approaches = list(results.keys())
    original = results[approaches[0]]
    improved = results[approaches[1]]

    print("\n" + "=" * 60)
    print("IMPROVEMENT ASSESSMENT")
    print("=" * 60)

    # Calculate improvements
    mean_improvement = (
        (original["mean_rf"] - improved["mean_rf"]) / original["mean_rf"]
    ) * 100
    max_improvement = (
        (original["max_rf"] - improved["max_rf"]) / original["max_rf"]
    ) * 100
    std_improvement = (
        (original["std_rf"] - improved["std_rf"]) / original["std_rf"]
    ) * 100

    print(f"Mean RF Distance:")
    print(f"  Original: {original['mean_rf']:.3f}")
    print(f"  Improved: {improved['mean_rf']:.3f}")
    print(f"  Improvement: {mean_improvement:+.1f}%")

    print(f"\nMaximum RF Distance:")
    print(f"  Original: {original['max_rf']:.3f}")
    print(f"  Improved: {improved['max_rf']:.3f}")
    print(f"  Improvement: {max_improvement:+.1f}%")

    print(f"\nStandard Deviation:")
    print(f"  Original: {original['std_rf']:.3f}")
    print(f"  Improved: {improved['std_rf']:.3f}")
    print(f"  Improvement: {std_improvement:+.1f}%")

    print(f"\nResolution:")
    print(f"  Original: {original['n_trees']} trees ({original['step_size']}bp steps)")
    print(f"  Improved: {improved['n_trees']} trees ({improved['step_size']}bp steps)")

    print(f"\nOverlap:")
    print(f"  Original: {original['overlap_pct']:.1f}%")
    print(f"  Improved: {improved['overlap_pct']:.1f}%")

    # Overall assessment
    print(f"\n🏆 OVERALL ASSESSMENT:")
    if mean_improvement > 10:
        print("✅ SIGNIFICANT IMPROVEMENT: Much better topological congruence")
    elif mean_improvement > 0:
        print("✅ MODERATE IMPROVEMENT: Better topological congruence")
    else:
        print("❌ NO IMPROVEMENT: Consider alternative strategies")

    if improved["mean_rf"] < 0.5:
        print("✅ EXCELLENT: Mean RF distance < 0.5 (very good congruence)")
    elif improved["mean_rf"] < 1.0:
        print("✅ GOOD: Mean RF distance < 1.0 (acceptable congruence)")
    else:
        print("⚠️  MODERATE: Still significant topological shifts")


def main():
    """Main comparison function."""

    # Compare approaches
    results = compare_approaches()

    if not results:
        print("No results to compare")
        return

    # Create comparison table
    create_comparison_table(results)

    # Create visualization
    create_comparison_visualization(results)

    # Assess improvement
    assess_improvement(results)


if __name__ == "__main__":
    main()
