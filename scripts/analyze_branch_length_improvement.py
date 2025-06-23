#!/usr/bin/env python3
"""
Analyze the branch lengths in the cleaned single outgroup results
to see if removing KR074191.1 and using a single outgroup resolved the issue.
"""

import os
import re
import statistics


def extract_branch_lengths(tree_string):
    """Extract all branch lengths from a Newick tree string."""
    # Find all numbers after colons (branch lengths)
    pattern = r":([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
    matches = re.findall(pattern, tree_string)
    return [float(match) for match in matches]


def analyze_tree_file(tree_file):
    """Analyze branch lengths in a tree file."""
    if not os.path.exists(tree_file):
        return None

    all_branch_lengths = []
    tree_stats = []

    with open(tree_file, "r") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if line:
                branch_lengths = extract_branch_lengths(line)
                if branch_lengths:
                    all_branch_lengths.extend(branch_lengths)

                    tree_stat = {
                        "tree_number": i,
                        "min_branch": min(branch_lengths),
                        "max_branch": max(branch_lengths),
                        "mean_branch": statistics.mean(branch_lengths),
                        "total_length": sum(branch_lengths),
                        "n_branches": len(branch_lengths),
                    }
                    tree_stats.append(tree_stat)

    if all_branch_lengths:
        overall_stats = {
            "total_branches": len(all_branch_lengths),
            "min_branch": min(all_branch_lengths),
            "max_branch": max(all_branch_lengths),
            "mean_branch": statistics.mean(all_branch_lengths),
            "median_branch": statistics.median(all_branch_lengths),
            "problematic_branches": sum(1 for bl in all_branch_lengths if bl > 5.0),
            "very_long_branches": sum(1 for bl in all_branch_lengths if bl > 1.0),
            "tree_stats": tree_stats,
        }
        return overall_stats

    return None


def main():
    """Compare branch lengths before and after cleanup."""
    print("=== BRANCH LENGTH ANALYSIS: BEFORE vs AFTER CLEANUP ===\n")

    # Define the result files to compare
    results_files = {
        "Original (with outgroups)": "/Users/berksakalli/Projects/automated-window-sliding/results_enhanced_outgroups/best_rooted_trees.newick",
        "Cleaned (single outgroup)": "/Users/berksakalli/Projects/automated-window-sliding/results_single_outgroup/best_rooted_trees.newick",
    }

    comparison_data = {}

    for run_name, file_path in results_files.items():
        print(f"--- Analyzing: {run_name} ---")

        if not os.path.exists(file_path):
            print(f"❌ File not found: {file_path}")
            continue

        stats = analyze_tree_file(file_path)

        if stats:
            comparison_data[run_name] = stats

            print(f"📊 Trees analyzed: {len(stats['tree_stats'])}")
            print(f"🌳 Total branches: {stats['total_branches']}")
            print(
                f"📏 Branch length range: {stats['min_branch']:.6f} - {stats['max_branch']:.6f}"
            )
            print(f"📈 Mean branch length: {stats['mean_branch']:.6f}")
            print(f"📊 Median branch length: {stats['median_branch']:.6f}")
            print(f"⚠️ Problematic branches (>5.0): {stats['problematic_branches']}")
            print(f"🔴 Very long branches (>1.0): {stats['very_long_branches']}")

            # Show worst trees
            worst_trees = sorted(
                stats["tree_stats"], key=lambda x: x["max_branch"], reverse=True
            )[:5]
            print(f"\n🔥 Top 5 trees with longest branches:")
            for tree in worst_trees:
                print(
                    f"   Tree {tree['tree_number']}: max={tree['max_branch']:.6f}, mean={tree['mean_branch']:.6f}"
                )
        else:
            print(f"❌ Could not analyze file: {file_path}")

        print()

    # Comparison summary
    if len(comparison_data) >= 2:
        original_key = "Original (with outgroups)"
        cleaned_key = "Cleaned (single outgroup)"

        if original_key in comparison_data and cleaned_key in comparison_data:
            print("=== IMPROVEMENT ANALYSIS ===")

            orig = comparison_data[original_key]
            clean = comparison_data[cleaned_key]

            print(
                f"📉 Max branch length reduction: {orig['max_branch']:.6f} → {clean['max_branch']:.6f}"
            )
            print(
                f"📉 Mean branch length change: {orig['mean_branch']:.6f} → {clean['mean_branch']:.6f}"
            )
            print(
                f"📉 Problematic branches (>5.0): {orig['problematic_branches']} → {clean['problematic_branches']}"
            )
            print(
                f"📉 Very long branches (>1.0): {orig['very_long_branches']} → {clean['very_long_branches']}"
            )

            # Calculate improvements
            max_improvement = (
                (orig["max_branch"] - clean["max_branch"]) / orig["max_branch"]
            ) * 100
            mean_improvement = (
                (orig["mean_branch"] - clean["mean_branch"]) / orig["mean_branch"]
            ) * 100
            problematic_reduction = (
                orig["problematic_branches"] - clean["problematic_branches"]
            )

            print(f"\n🎯 IMPROVEMENTS:")
            print(f"   Max branch length reduced by {max_improvement:.1f}%")
            print(f"   Mean branch length reduced by {mean_improvement:.1f}%")
            print(f"   Removed {problematic_reduction} problematic branches")

            if clean["max_branch"] < 1.0:
                print(
                    "✅ SUCCESS: Maximum branch length now under 1.0 (reasonable range)"
                )
            elif clean["max_branch"] < 2.0:
                print("✅ GOOD: Maximum branch length under 2.0 (acceptable range)")
            else:
                print(
                    "⚠️ CAUTION: Maximum branch length still high, may need further optimization"
                )

    # Check sequence count
    print("\n=== SEQUENCE COUNT VERIFICATION ===")

    # Check the alignment used
    single_outgroup_alignment = "/Users/berksakalli/Projects/automated-window-sliding/single_outgroup_alignment.fasta"
    if os.path.exists(single_outgroup_alignment):
        seq_count = 0
        with open(single_outgroup_alignment, "r") as f:
            for line in f:
                if line.startswith(">"):
                    seq_count += 1
        print(f"✅ Single outgroup alignment: {seq_count} sequences")

        # Check if KR074191.1 was removed
        with open(single_outgroup_alignment, "r") as f:
            content = f.read()
            if "KR074191" in content:
                print("⚠️ WARNING: KR074191.1 still present in alignment")
            else:
                print("✅ CONFIRMED: KR074191.1 successfully removed")

            # Check outgroup count
            outgroup_count = sum(
                1 for og in ["M87661", "L07418", "AF093797"] if og in content
            )
            print(
                f"🎯 Outgroups present: {outgroup_count}/3 (expecting 1 for single outgroup strategy)"
            )
    else:
        print(f"❌ Alignment file not found: {single_outgroup_alignment}")


if __name__ == "__main__":
    main()
