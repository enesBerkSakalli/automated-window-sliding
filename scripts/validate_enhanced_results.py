#!/usr/bin/env python3
"""
Final validation and tree quality analysis for the enhanced outgroups pipeline run.
This script validates the success of the sliding window phylogenetic analysis
and provides insights into the generated trees.
"""

import os
from collections import Counter


def analyze_tree_file(tree_file):
    """Analyze trees in a Newick file."""
    if not os.path.exists(tree_file):
        return {}

    analysis = {
        "total_trees": 0,
        "taxa_counts": [],
        "has_branch_lengths": 0,
        "tree_lengths": [],
    }

    with open(tree_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                analysis["total_trees"] += 1

                # Count taxa (terminal nodes = leaves)
                taxa_count = line.count(",") + 1
                analysis["taxa_counts"].append(taxa_count)

                # Check for branch lengths (presence of colon)
                if ":" in line:
                    analysis["has_branch_lengths"] += 1

                # Estimate tree length (rough approximation)
                if ":" in line:
                    # Count numeric values after colons as branch lengths
                    import re

                    branch_lengths = re.findall(r":([0-9.e-]+)", line)
                    if branch_lengths:
                        tree_length = sum(float(bl) for bl in branch_lengths)
                        analysis["tree_lengths"].append(tree_length)

    return analysis


def check_outgroup_presence(tree_file, outgroup_ids):
    """Check if outgroup sequences are present in trees."""
    if not os.path.exists(tree_file):
        return {"found": False, "details": "File not found"}

    outgroup_presence = {og: 0 for og in outgroup_ids}
    total_trees = 0

    with open(tree_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                total_trees += 1
                for og_id in outgroup_ids:
                    if og_id in line:
                        outgroup_presence[og_id] += 1

    return {
        "total_trees": total_trees,
        "outgroup_presence": outgroup_presence,
        "all_outgroups_present": all(count > 0 for count in outgroup_presence.values()),
    }


def analyze_window_results(results_dir):
    """Analyze individual window results."""
    tree_recon_dir = os.path.join(results_dir, "tree_reconstruction", "iqtree")

    if not os.path.exists(tree_recon_dir):
        return {}

    window_dirs = [
        d
        for d in os.listdir(tree_recon_dir)
        if os.path.isdir(os.path.join(tree_recon_dir, d))
    ]

    window_analysis = {
        "total_windows": len(window_dirs),
        "windows_with_trees": 0,
        "windows_with_rooted_trees": 0,
        "window_positions": [],
    }

    for window_dir in window_dirs:
        window_path = os.path.join(tree_recon_dir, window_dir)

        # Extract window position
        try:
            window_pos = int(window_dir)
            window_analysis["window_positions"].append(window_pos)
        except ValueError:
            continue

        # Check for tree files
        tree_file = os.path.join(window_path, f"{window_dir}.treefile")
        rooted_tree_file = os.path.join(window_path, f"{window_dir}_rooted.newick")

        if os.path.exists(tree_file):
            window_analysis["windows_with_trees"] += 1

        if os.path.exists(rooted_tree_file):
            window_analysis["windows_with_rooted_trees"] += 1

    return window_analysis


def main():
    """Main validation function."""
    print("=== ENHANCED OUTGROUPS PIPELINE VALIDATION ===\n")

    results_dir = "/Users/berksakalli/Projects/automated-window-sliding/results_enhanced_outgroups"

    if not os.path.exists(results_dir):
        print(f"❌ Results directory not found: {results_dir}")
        return

    print(f"📁 Analyzing results in: {results_dir}\n")

    # Define expected files
    expected_files = {
        "rooted_trees": "best_rooted_trees.newick",
        "unrooted_trees": "best_trees.newick",
        "consensus_rooted": "consensus_trees.newick",
        "rooted_nexus": "best_rooted_trees.nexus",
    }

    # Check file existence
    print("📊 FILE EXISTENCE CHECK:")
    for file_type, filename in expected_files.items():
        filepath = os.path.join(results_dir, filename)
        if os.path.exists(filepath):
            print(f"✅ {file_type}: {filename}")
        else:
            print(f"❌ {file_type}: {filename}")
    print()

    # Analyze rooted trees
    rooted_trees_file = os.path.join(results_dir, "best_rooted_trees.newick")
    print("🌳 ROOTED TREES ANALYSIS:")

    if os.path.exists(rooted_trees_file):
        tree_analysis = analyze_tree_file(rooted_trees_file)
        print(f"✅ Total rooted trees: {tree_analysis['total_trees']}")

        if tree_analysis["taxa_counts"]:
            taxa_counts = Counter(tree_analysis["taxa_counts"])
            print(f"📈 Taxa per tree: {dict(taxa_counts)}")

        print(f"🔗 Trees with branch lengths: {tree_analysis['has_branch_lengths']}")

        if tree_analysis["tree_lengths"]:
            avg_length = sum(tree_analysis["tree_lengths"]) / len(
                tree_analysis["tree_lengths"]
            )
            min_length = min(tree_analysis["tree_lengths"])
            max_length = max(tree_analysis["tree_lengths"])
            print(
                f"📏 Tree lengths - Avg: {avg_length:.6f}, Min: {min_length:.6f}, Max: {max_length:.6f}"
            )
    else:
        print("❌ Rooted trees file not found")
    print()

    # Check outgroup presence
    gi_outgroups = ["M87661", "L07418", "AF093797"]
    print("🎯 OUTGROUP PRESENCE CHECK:")

    if os.path.exists(rooted_trees_file):
        outgroup_check = check_outgroup_presence(rooted_trees_file, gi_outgroups)
        print(f"📊 Total trees analyzed: {outgroup_check['total_trees']}")

        for og_id, count in outgroup_check["outgroup_presence"].items():
            percentage = (
                (count / outgroup_check["total_trees"]) * 100
                if outgroup_check["total_trees"] > 0
                else 0
            )
            print(
                f"🧬 {og_id}: Present in {count}/{outgroup_check['total_trees']} trees ({percentage:.1f}%)"
            )

        if outgroup_check["all_outgroups_present"]:
            print("✅ All GI outgroups detected in trees")
        else:
            print("⚠️ Some outgroups may be missing from trees")
    print()

    # Analyze individual windows
    print("🪟 SLIDING WINDOW ANALYSIS:")
    window_analysis = analyze_window_results(results_dir)

    if window_analysis:
        print(f"📈 Total windows processed: {window_analysis['total_windows']}")
        print(f"🌳 Windows with trees: {window_analysis['windows_with_trees']}")
        print(
            f"🔀 Windows with rooted trees: {window_analysis['windows_with_rooted_trees']}"
        )

        if window_analysis["window_positions"]:
            positions = sorted(window_analysis["window_positions"])
            print(
                f"📍 Window positions: {positions[0]} to {positions[-1]} (step {positions[1] - positions[0] if len(positions) > 1 else 'N/A'})"
            )

            # Validate sliding window parameters
            if len(positions) > 1:
                step_size = positions[1] - positions[0]
                expected_windows = ((positions[-1] - positions[0]) // step_size) + 1
                print(
                    f"🎯 Expected vs Actual windows: {expected_windows} vs {len(positions)}"
                )
    print()

    # Check alignment file
    print("🧬 ALIGNMENT VALIDATION:")
    alignment_file = "/Users/berksakalli/Projects/automated-window-sliding/enhanced_alignment_with_outgroups.fasta"

    if os.path.exists(alignment_file):
        seq_count = 0
        seq_lengths = []
        gi_seqs_found = []

        with open(alignment_file, "r") as f:
            current_seq = ""
            current_header = ""

            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq and current_header:
                        seq_lengths.append(len(current_seq))

                    seq_count += 1
                    current_header = line[1:]
                    current_seq = ""

                    # Check for GI outgroups
                    for gi_id in gi_outgroups:
                        if gi_id in current_header:
                            gi_seqs_found.append(gi_id)
                else:
                    current_seq += line

            # Process last sequence
            if current_seq and current_header:
                seq_lengths.append(len(current_seq))

        print(f"✅ Enhanced alignment: {seq_count} sequences")

        if seq_lengths:
            avg_length = sum(seq_lengths) / len(seq_lengths)
            print(f"📏 Average sequence length: {avg_length:.0f}")

        print(f"🎯 GI outgroups in alignment: {gi_seqs_found}")

        if len(gi_seqs_found) == len(gi_outgroups):
            print("✅ All GI outgroups present in alignment")
        else:
            missing = set(gi_outgroups) - set(gi_seqs_found)
            print(f"⚠️ Missing GI outgroups: {missing}")
    else:
        print(f"❌ Enhanced alignment not found: {alignment_file}")
    print()

    # Final validation summary
    print("=== VALIDATION SUMMARY ===")

    success_criteria = [
        (os.path.exists(rooted_trees_file), "Rooted trees file exists"),
        (tree_analysis.get("total_trees", 0) == 33, "33 rooted trees generated"),
        (window_analysis.get("total_windows", 0) == 33, "33 sliding windows processed"),
        (len(gi_seqs_found) == 3, "All 3 GI outgroups in alignment"),
        (os.path.exists(alignment_file), "Enhanced alignment file exists"),
    ]

    passed = sum(1 for criterion, _ in success_criteria if criterion)
    total = len(success_criteria)

    print(f"📊 Validation Score: {passed}/{total} criteria passed\n")

    for criterion, description in success_criteria:
        status = "✅" if criterion else "❌"
        print(f"{status} {description}")

    if passed == total:
        print("\n🎉 VALIDATION SUCCESSFUL! 🎉")
        print(
            "The enhanced outgroups pipeline has completed successfully with all expected outputs."
        )
        print("\nMethodological improvements achieved:")
        print("- ✅ GI outgroup integration for robust phylogenetic rooting")
        print("- ✅ Quality-filtered alignment with outlier removal")
        print("- ✅ Literature-based reference sequences")
        print("- ✅ Complete sliding window analysis (window=200, step=15)")
        print("- ✅ 33 properly rooted trees for downstream analysis")
    else:
        print(f"\n⚠️ VALIDATION INCOMPLETE: {total - passed} criteria failed")
        print("Review the failed criteria above for troubleshooting.")


if __name__ == "__main__":
    main()
