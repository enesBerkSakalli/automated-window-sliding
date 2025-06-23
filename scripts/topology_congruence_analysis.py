#!/usr/bin/env python3
"""
Tree Topology Congruence Analysis
Analyzes topological differences between adjacent sliding windows and suggests improvements.
"""

import os
import numpy as np
import dendropy
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt
import pandas as pd


def load_trees_from_newick(file_path):
    """Load all trees from a Newick file."""
    trees = []

    with open(file_path, "r") as f:
        content = f.read().strip()

    tree_lines = [line.strip() for line in content.split("\n") if line.strip()]

    for i, tree_line in enumerate(tree_lines):
        if tree_line:
            try:
                # Use DendroPy for more robust tree comparison
                tree = dendropy.Tree.get(data=tree_line, schema="newick")
                trees.append(tree)
                print(f"Loaded tree {i + 1}: {len(tree.leaf_nodes())} taxa")
            except Exception as e:
                print(f"Error loading tree {i + 1}: {e}")

    return trees


def calculate_robinson_foulds_distances(trees):
    """Calculate Robinson-Foulds distances between all pairs of trees."""
    n_trees = len(trees)
    rf_matrix = np.zeros((n_trees, n_trees))

    print(f"Calculating Robinson-Foulds distances for {n_trees} trees...")

    # Ensure all trees have the same taxon namespace
    taxon_namespace = trees[0].taxon_namespace
    for tree in trees[1:]:
        tree.migrate_taxon_namespace(taxon_namespace)

    for i in range(n_trees):
        for j in range(i, n_trees):
            if i == j:
                rf_matrix[i, j] = 0
            else:
                try:
                    rf_distance = treecompare.robinson_foulds_distance(
                        trees[i], trees[j]
                    )
                    rf_matrix[i, j] = rf_distance
                    rf_matrix[j, i] = rf_distance
                except Exception as e:
                    print(
                        f"Error calculating RF distance between trees {i + 1} and {j + 1}: {e}"
                    )
                    rf_matrix[i, j] = np.nan
                    rf_matrix[j, i] = np.nan

    return rf_matrix


def analyze_adjacent_differences(rf_matrix):
    """Analyze differences between adjacent sliding windows."""
    n_trees = len(rf_matrix)
    adjacent_distances = []

    for i in range(n_trees - 1):
        distance = rf_matrix[i, i + 1]
        if not np.isnan(distance):
            adjacent_distances.append(distance)
            print(f"Window {i + 1} → {i + 2}: RF distance = {distance}")

    if adjacent_distances:
        print(f"\nAdjacent window RF distances:")
        print(f"Mean: {np.mean(adjacent_distances):.2f}")
        print(f"Median: {np.median(adjacent_distances):.2f}")
        print(f"Min: {np.min(adjacent_distances):.2f}")
        print(f"Max: {np.max(adjacent_distances):.2f}")
        print(f"Std: {np.std(adjacent_distances):.2f}")

    return adjacent_distances


def identify_problematic_transitions(adjacent_distances, threshold_percentile=75):
    """Identify windows with unusually large topological changes."""
    if not adjacent_distances:
        return []

    threshold = np.percentile(adjacent_distances, threshold_percentile)
    problematic = []

    print(f"\nProblematic transitions (>{threshold:.1f} RF distance):")
    for i, distance in enumerate(adjacent_distances):
        if distance > threshold:
            problematic.append((i + 1, i + 2, distance))
            print(f"Window {i + 1} → {i + 2}: RF = {distance:.1f} (HIGH)")

    return problematic


def visualize_topology_changes(
    rf_matrix, adjacent_distances, output_file="topology_analysis.png"
):
    """Create visualizations of topological changes."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    # 1. RF distance matrix heatmap
    im1 = axes[0, 0].imshow(rf_matrix, cmap="viridis", aspect="auto")
    axes[0, 0].set_title("Robinson-Foulds Distance Matrix")
    axes[0, 0].set_xlabel("Tree Number")
    axes[0, 0].set_ylabel("Tree Number")
    plt.colorbar(im1, ax=axes[0, 0])

    # 2. Adjacent window distances
    if adjacent_distances:
        axes[0, 1].plot(
            range(1, len(adjacent_distances) + 1), adjacent_distances, "o-", color="red"
        )
        axes[0, 1].set_title("Adjacent Window RF Distances")
        axes[0, 1].set_xlabel("Window Transition")
        axes[0, 1].set_ylabel("RF Distance")
        axes[0, 1].axhline(
            np.mean(adjacent_distances),
            color="blue",
            linestyle="--",
            label=f"Mean: {np.mean(adjacent_distances):.1f}",
        )
        axes[0, 1].legend()

    # 3. Distribution of adjacent distances
    if adjacent_distances:
        axes[1, 0].hist(
            adjacent_distances, bins=10, alpha=0.7, color="skyblue", edgecolor="black"
        )
        axes[1, 0].set_title("Distribution of Adjacent RF Distances")
        axes[1, 0].set_xlabel("RF Distance")
        axes[1, 0].set_ylabel("Frequency")
        axes[1, 0].axvline(
            np.mean(adjacent_distances),
            color="red",
            linestyle="--",
            label=f"Mean: {np.mean(adjacent_distances):.1f}",
        )
        axes[1, 0].legend()

    # 4. Cumulative distribution
    if adjacent_distances:
        sorted_distances = np.sort(adjacent_distances)
        cumulative = np.arange(1, len(sorted_distances) + 1) / len(sorted_distances)
        axes[1, 1].plot(sorted_distances, cumulative, color="purple", linewidth=2)
        axes[1, 1].set_title("Cumulative Distribution of RF Distances")
        axes[1, 1].set_xlabel("RF Distance")
        axes[1, 1].set_ylabel("Cumulative Probability")
        axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Topology analysis saved as {output_file}")


def suggest_improvements():
    """Suggest strategies to improve topological congruence."""
    print("\n" + "=" * 70)
    print("STRATEGIES TO IMPROVE TOPOLOGICAL CONGRUENCE")
    print("=" * 70)

    print("\n1. 🔧 PARAMETER ADJUSTMENTS:")
    print("-" * 35)
    print("   • Increase window size (200bp → 300-500bp)")
    print("     - Larger windows = more phylogenetic signal")
    print("     - Reduces noise from alignment artifacts")
    print("   • Decrease step size (15bp → 5-10bp)")
    print("     - Smaller steps = smoother transitions")
    print("     - More gradual topological changes")
    print("   • Increase minimum window overlap")
    print("     - Ensures continuity between adjacent windows")

    print("\n2. 🧬 ALIGNMENT IMPROVEMENTS:")
    print("-" * 36)
    print("   • Use codon-aware alignment (MACSE, PRANK)")
    print("     - Maintains reading frame integrity")
    print("     - Reduces alignment artifacts")
    print("   • Apply stricter quality filtering")
    print("     - Remove highly variable regions")
    print("     - Focus on conserved phylogenetic markers")
    print("   • Use profile-based alignment")
    print("     - More consistent gap placement")
    print("     - Better handling of indels")

    print("\n3. 🌳 PHYLOGENETIC CONSTRAINTS:")
    print("-" * 38)
    print("   • Apply topological constraints")
    print("     - Use backbone tree from full genome")
    print("     - Constrain major clades")
    print("   • Use consensus tree as guide")
    print("     - Start each window with similar topology")
    print("     - Reduce search space")
    print("   • Implement tree smoothing")
    print("     - Weight trees by adjacent window similarity")
    print("     - Use sliding window consensus")

    print("\n4. 📊 STATISTICAL APPROACHES:")
    print("-" * 37)
    print("   • Use bootstrap consensus trees")
    print("     - More stable topologies")
    print("     - Reduced noise from weak branches")
    print("   • Apply branch length smoothing")
    print("     - Gradual changes between windows")
    print("     - Reduce phylogenetic artifacts")
    print("   • Implement Bayesian sliding windows")
    print("     - Posterior probability weighting")
    print("     - More conservative topology changes")

    print("\n5. 🔄 POST-PROCESSING:")
    print("-" * 26)
    print("   • Tree reconciliation")
    print("     - Minimize total RF distance")
    print("     - Optimal tree ordering")
    print("   • Consensus approaches")
    print("     - Extended majority rule consensus")
    print("     - Weighted by window overlap")
    print("   • Outlier detection and smoothing")
    print("     - Identify aberrant topologies")
    print("     - Replace with interpolated trees")


def analyze_sequence_properties():
    """Analyze properties that might cause topological instability."""
    print("\n" + "=" * 70)
    print("SEQUENCE PROPERTIES ANALYSIS")
    print("=" * 70)

    alignment_file = "ingroup_only_alignment.fasta"
    if not os.path.exists(alignment_file):
        print(f"Alignment file {alignment_file} not found")
        return

    from Bio import SeqIO, AlignIO

    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        seq_length = alignment.get_alignment_length()
        n_sequences = len(alignment)

        print(f"Alignment properties:")
        print(f"• Sequences: {n_sequences}")
        print(f"• Length: {seq_length} bp")
        print(f"• Window size: 200bp ({200 / seq_length * 100:.1f}% of genome)")
        print(f"• Number of windows: 33")
        print(f"• Step size: 15bp")
        print(f"• Window overlap: {(200 - 15) / 200 * 100:.1f}%")

        # Calculate variability across the alignment
        variability = []
        for i in range(seq_length):
            column = alignment[:, i]
            unique_chars = len(set(column.replace("-", "").replace("N", "")))
            variability.append(unique_chars)

        # Analyze variability in windows
        window_variability = []
        for start in range(0, seq_length - 200 + 1, 15):
            end = start + 200
            window_var = np.mean(variability[start:end])
            window_variability.append(window_var)

        if window_variability:
            print(f"\nWindow variability analysis:")
            print(f"• Mean variability per window: {np.mean(window_variability):.2f}")
            print(
                f"• Variability range: {np.min(window_variability):.2f} - {np.max(window_variability):.2f}"
            )
            print(
                f"• High variability windows: {np.sum(np.array(window_variability) > np.percentile(window_variability, 75))}"
            )

    except Exception as e:
        print(f"Error analyzing sequence properties: {e}")


def main():
    """Main analysis function."""
    print("Tree Topology Congruence Analysis")
    print("=" * 50)

    # Load trees
    trees_file = "results_midpoint_rooting/best_rooted_trees.newick"
    if not os.path.exists(trees_file):
        print(f"Trees file {trees_file} not found")
        return

    print(f"Loading trees from {trees_file}...")
    trees = load_trees_from_newick(trees_file)

    if len(trees) < 2:
        print("Need at least 2 trees for comparison")
        return

    # Calculate Robinson-Foulds distances
    rf_matrix = calculate_robinson_foulds_distances(trees)

    # Analyze adjacent differences
    adjacent_distances = analyze_adjacent_differences(rf_matrix)

    # Identify problematic transitions
    problematic = identify_problematic_transitions(adjacent_distances)

    # Create visualizations
    visualize_topology_changes(rf_matrix, adjacent_distances)

    # Analyze sequence properties
    analyze_sequence_properties()

    # Suggest improvements
    suggest_improvements()

    # Save results
    if adjacent_distances:
        df = pd.DataFrame(
            {
                "Window_Transition": [
                    f"{i + 1}→{i + 2}" for i in range(len(adjacent_distances))
                ],
                "RF_Distance": adjacent_distances,
            }
        )
        df.to_csv("topology_congruence_analysis.csv", index=False)
        print(f"\nResults saved to topology_congruence_analysis.csv")


if __name__ == "__main__":
    main()
