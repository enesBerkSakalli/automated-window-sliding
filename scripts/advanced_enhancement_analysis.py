#!/usr/bin/env python3
"""
Advanced Topological Congruence Analysis
Deep analysis of current results and identification of further improvements.
"""

import os
import numpy as np
import dendropy
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO
import seaborn as sns


def detailed_rf_analysis(trees_file, approach_name):
    """Perform detailed Robinson-Foulds analysis."""

    print(f"\n{'=' * 60}")
    print(f"DETAILED ANALYSIS: {approach_name}")
    print(f"{'=' * 60}")

    if not os.path.exists(trees_file):
        print(f"File not found: {trees_file}")
        return None

    # Load trees
    trees = []
    with open(trees_file, "r") as f:
        content = f.read().strip()

    tree_lines = [line.strip() for line in content.split("\n") if line.strip()]

    for tree_line in tree_lines:
        if tree_line:
            tree = dendropy.Tree.get(data=tree_line, schema="newick")
            trees.append(tree)

    n_trees = len(trees)
    print(f"Loaded {n_trees} trees")

    # Ensure same taxon namespace
    taxon_namespace = trees[0].taxon_namespace
    for tree in trees[1:]:
        tree.migrate_taxon_namespace(taxon_namespace)

    # Calculate all pairwise RF distances
    rf_matrix = np.zeros((n_trees, n_trees))

    for i in range(n_trees):
        for j in range(i + 1, n_trees):
            rf_dist = treecompare.robinson_foulds_distance(trees[i], trees[j])
            rf_matrix[i, j] = rf_dist
            rf_matrix[j, i] = rf_dist

    # Adjacent window analysis
    adjacent_distances = []
    for i in range(n_trees - 1):
        adjacent_distances.append(rf_matrix[i, i + 1])

    # Sliding window analysis (3-window, 5-window, 10-window)
    analyses = {}

    for window_size in [3, 5, 10]:
        if window_size <= n_trees:
            window_distances = []
            for i in range(n_trees - window_size + 1):
                window_trees = trees[i : i + window_size]

                # Calculate mean RF within window
                window_rf = []
                for j in range(len(window_trees)):
                    for k in range(j + 1, len(window_trees)):
                        rf = treecompare.robinson_foulds_distance(
                            window_trees[j], window_trees[k]
                        )
                        window_rf.append(rf)

                if window_rf:
                    window_distances.append(np.mean(window_rf))

            analyses[f"{window_size}-window"] = window_distances

    # Find problematic regions
    threshold_75 = np.percentile(adjacent_distances, 75)
    threshold_90 = np.percentile(adjacent_distances, 90)

    problematic_regions = []
    for i, rf_dist in enumerate(adjacent_distances):
        if rf_dist > threshold_90:
            problematic_regions.append(
                {
                    "transition": f"{i + 1}→{i + 2}",
                    "rf_distance": rf_dist,
                    "severity": "High",
                }
            )
        elif rf_dist > threshold_75:
            problematic_regions.append(
                {
                    "transition": f"{i + 1}→{i + 2}",
                    "rf_distance": rf_dist,
                    "severity": "Moderate",
                }
            )

    # Statistics
    stats = {
        "n_trees": n_trees,
        "adjacent_mean": np.mean(adjacent_distances),
        "adjacent_median": np.median(adjacent_distances),
        "adjacent_std": np.std(adjacent_distances),
        "adjacent_min": np.min(adjacent_distances),
        "adjacent_max": np.max(adjacent_distances),
        "problematic_high": len(
            [p for p in problematic_regions if p["severity"] == "High"]
        ),
        "problematic_moderate": len(
            [p for p in problematic_regions if p["severity"] == "Moderate"]
        ),
        "rf_matrix": rf_matrix,
        "adjacent_distances": adjacent_distances,
        "window_analyses": analyses,
        "problematic_regions": problematic_regions,
    }

    print(f"\nAdjacent Window Statistics:")
    print(f"  Mean RF: {stats['adjacent_mean']:.3f}")
    print(f"  Median RF: {stats['adjacent_median']:.3f}")
    print(f"  Std RF: {stats['adjacent_std']:.3f}")
    print(f"  Range: {stats['adjacent_min']:.3f} - {stats['adjacent_max']:.3f}")
    print(f"  High problems (>90th percentile): {stats['problematic_high']}")
    print(f"  Moderate problems (>75th percentile): {stats['problematic_moderate']}")

    return stats


def analyze_rooting_consistency(trees_file, approach_name):
    """Analyze rooting consistency across trees."""

    print(f"\n{'=' * 60}")
    print(f"ROOTING CONSISTENCY ANALYSIS: {approach_name}")
    print(f"{'=' * 60}")

    if not os.path.exists(trees_file):
        print(f"File not found: {trees_file}")
        return None

    trees = []
    with open(trees_file, "r") as f:
        content = f.read().strip()

    tree_lines = [line.strip() for line in content.split("\n") if line.strip()]

    root_positions = []
    tree_stats = []

    for i, tree_line in enumerate(tree_lines):
        if tree_line:
            tree = dendropy.Tree.get(data=tree_line, schema="newick")
            trees.append(tree)

            # Analyze root position
            if tree.seed_node and len(tree.seed_node.child_nodes()) >= 2:
                left_child = tree.seed_node.child_nodes()[0]
                right_child = tree.seed_node.child_nodes()[1]

                left_taxa = len(left_child.leaf_nodes())
                right_taxa = len(right_child.leaf_nodes())
                total_taxa = left_taxa + right_taxa

                # Calculate balance (closer to 0.5 = more balanced)
                balance = min(left_taxa, right_taxa) / total_taxa
                root_positions.append(balance)

                # Tree statistics
                tree_length = tree.length()
                max_dist = tree.max_distance_from_root()

                tree_stats.append(
                    {
                        "tree_id": i + 1,
                        "balance": balance,
                        "tree_length": tree_length,
                        "max_distance": max_dist,
                        "left_taxa": left_taxa,
                        "right_taxa": right_taxa,
                    }
                )

    if root_positions:
        print(f"Root Balance Analysis ({len(root_positions)} trees):")
        print(f"  Mean balance: {np.mean(root_positions):.3f}")
        print(f"  Std balance: {np.std(root_positions):.3f}")
        print(f"  Range: {np.min(root_positions):.3f} - {np.max(root_positions):.3f}")

        # Check for consistency
        balance_std = np.std(root_positions)
        if balance_std < 0.05:
            print("  ✅ Excellent rooting consistency")
        elif balance_std < 0.1:
            print("  ✅ Good rooting consistency")
        elif balance_std < 0.2:
            print("  ⚠️  Moderate rooting consistency")
        else:
            print("  ❌ Poor rooting consistency")

    return {
        "root_positions": root_positions,
        "tree_stats": tree_stats,
        "balance_mean": np.mean(root_positions) if root_positions else 0,
        "balance_std": np.std(root_positions) if root_positions else 0,
    }


def identify_improvement_opportunities():
    """Identify specific opportunities for further improvement."""

    print(f"\n{'=' * 60}")
    print("IMPROVEMENT OPPORTUNITIES ANALYSIS")
    print(f"{'=' * 60}")

    # Analyze current results
    current_results = detailed_rf_analysis(
        "results_improved_strategy2/best_rooted_trees.newick", "Current Improved"
    )

    rooting_analysis = analyze_rooting_consistency(
        "results_improved_strategy2/best_rooted_trees.newick", "Current Improved"
    )

    if not current_results:
        print("Cannot analyze current results")
        return

    improvements = []

    # 1. RF Distance Analysis
    if current_results["adjacent_mean"] > 0.5:
        improvements.append(
            {
                "category": "Topological Stability",
                "issue": f"Mean RF distance still high: {current_results['adjacent_mean']:.3f}",
                "solutions": [
                    "Further increase window size (400bp → 500-600bp)",
                    "Reduce step size even more (5bp → 3bp)",
                    "Apply topological constraints",
                    "Use consensus tree approaches",
                ],
                "priority": "High"
                if current_results["adjacent_mean"] > 1.0
                else "Medium",
            }
        )

    # 2. Problematic Regions
    if current_results["problematic_high"] > 0:
        improvements.append(
            {
                "category": "Problematic Transitions",
                "issue": f"{current_results['problematic_high']} high-problem transitions",
                "solutions": [
                    "Apply outlier tree smoothing",
                    "Use weighted consensus for problematic regions",
                    "Implement tree reconciliation algorithms",
                    "Apply bootstrap consensus thresholds",
                ],
                "priority": "High",
            }
        )

    # 3. Rooting Consistency
    if rooting_analysis and rooting_analysis["balance_std"] > 0.1:
        improvements.append(
            {
                "category": "Rooting Consistency",
                "issue": f"Variable rooting balance (std: {rooting_analysis['balance_std']:.3f})",
                "solutions": [
                    "Use constrained rooting with backbone tree",
                    "Apply multiple outgroup approach",
                    "Implement midpoint rooting consistency checks",
                    "Use Bayesian rooting with priors",
                ],
                "priority": "High",
            }
        )

    # 4. Resolution vs Quality Trade-off
    improvements.append(
        {
            "category": "Resolution Enhancement",
            "issue": "Balance between resolution and quality",
            "solutions": [
                "Test 2bp steps for maximum resolution",
                "Implement adaptive step sizes",
                "Use overlapping consensus windows",
                "Apply spline smoothing to tree space",
            ],
            "priority": "Medium",
        }
    )

    # 5. Statistical Robustness
    improvements.append(
        {
            "category": "Statistical Robustness",
            "issue": "Single tree per window may be unstable",
            "solutions": [
                "Generate bootstrap consensus trees",
                "Use Bayesian posterior tree sets",
                "Apply model averaging approaches",
                "Implement uncertainty quantification",
            ],
            "priority": "Medium",
        }
    )

    return improvements


def create_enhancement_strategies():
    """Create specific enhancement strategies."""

    print(f"\n{'=' * 60}")
    print("ENHANCEMENT STRATEGIES")
    print(f"{'=' * 60}")

    strategies = {
        "Ultra-Fine Resolution": {
            "window_size": 450,
            "step_size": 2,
            "overlap_pct": 99.6,
            "description": "Maximum resolution with 2bp steps",
            "pros": [
                "Highest resolution",
                "Smoothest transitions",
                "Best for fine recombination mapping",
            ],
            "cons": ["More computationally intensive", "Many more trees to analyze"],
            "recommended_for": "Publication-quality fine-scale analysis",
        },
        "Constrained Topology": {
            "window_size": 400,
            "step_size": 5,
            "overlap_pct": 98.8,
            "use_constraints": True,
            "description": "Current parameters + topological constraints",
            "pros": [
                "Maintains current resolution",
                "Improved topology stability",
                "Faster convergence",
            ],
            "cons": ["May miss some recombination events", "Requires backbone tree"],
            "recommended_for": "Balanced approach with constraints",
        },
        "Bootstrap Consensus": {
            "window_size": 400,
            "step_size": 5,
            "overlap_pct": 98.8,
            "bootstrap_replicates": 100,
            "description": "Bootstrap consensus trees per window",
            "pros": ["Statistical robustness", "Confidence measures", "Reduced noise"],
            "cons": ["Much more computation", "Conservative results"],
            "recommended_for": "High-confidence publication results",
        },
        "Adaptive Windows": {
            "description": "Variable window sizes based on sequence variability",
            "pros": [
                "Optimized for local sequence properties",
                "Better signal in variable regions",
            ],
            "cons": ["Complex implementation", "Non-uniform analysis"],
            "recommended_for": "Advanced research applications",
        },
        "Multi-Scale Analysis": {
            "description": "Multiple window sizes analyzed in parallel",
            "window_sizes": [300, 400, 500],
            "step_size": 5,
            "pros": [
                "Multiple resolution levels",
                "Cross-validation",
                "Robust detection",
            ],
            "cons": ["3x computation", "Complex result interpretation"],
            "recommended_for": "Comprehensive recombination analysis",
        },
    }

    return strategies


def main():
    """Main enhancement analysis."""

    print("Advanced Topological Congruence Enhancement Analysis")
    print("=" * 60)

    # Detailed analysis of current results
    current_stats = detailed_rf_analysis(
        "results_improved_strategy2/best_rooted_trees.newick",
        "Current Improved (400bp/5bp)",
    )

    # Rooting analysis
    rooting_stats = analyze_rooting_consistency(
        "results_improved_strategy2/best_rooted_trees.newick", "Current Improved"
    )

    # Compare with original
    original_stats = detailed_rf_analysis(
        "results_midpoint_rooting/best_rooted_trees.newick", "Original (200bp/15bp)"
    )

    # Identify improvement opportunities
    improvements = identify_improvement_opportunities()

    print(f"\n{'=' * 60}")
    print("SPECIFIC IMPROVEMENT RECOMMENDATIONS")
    print(f"{'=' * 60}")

    for i, improvement in enumerate(improvements, 1):
        print(f"\n{i}. {improvement['category']} [{improvement['priority']} Priority]")
        print(f"   Issue: {improvement['issue']}")
        print(f"   Solutions:")
        for solution in improvement["solutions"]:
            print(f"     • {solution}")

    # Enhancement strategies
    strategies = create_enhancement_strategies()

    print(f"\n{'=' * 60}")
    print("RECOMMENDED ENHANCEMENT STRATEGIES")
    print(f"{'=' * 60}")

    for name, strategy in strategies.items():
        print(f"\n🔧 {name}:")
        print(f"   {strategy['description']}")
        if "window_size" in strategy:
            print(
                f"   Parameters: {strategy['window_size']}bp windows, {strategy['step_size']}bp steps"
            )
            print(f"   Overlap: {strategy['overlap_pct']:.1f}%")
        print(f"   ✅ Pros: {', '.join(strategy['pros'])}")
        print(f"   ⚠️  Cons: {', '.join(strategy['cons'])}")
        print(f"   🎯 Best for: {strategy['recommended_for']}")


if __name__ == "__main__":
    main()
