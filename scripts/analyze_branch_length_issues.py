#!/usr/bin/env python3
"""
Analyze and fix high branch length issues in phylogenetic trees.
This script identifies problematic sequences and provides solutions.
"""

import os
import re
from collections import defaultdict
import numpy as np


def extract_branch_lengths_from_tree(tree_string):
    """Extract all branch lengths from a Newick tree string."""
    # Find all branch lengths (numbers after colons)
    branch_lengths = re.findall(r":([0-9.e-]+)", tree_string)
    return [float(bl) for bl in branch_lengths]


def analyze_tree_branch_lengths(tree_file):
    """Analyze branch lengths across all trees."""
    if not os.path.exists(tree_file):
        return {}

    all_branch_lengths = []
    tree_stats = []
    problematic_trees = []

    with open(tree_file, "r") as f:
        for tree_idx, line in enumerate(f, 1):
            line = line.strip()
            if line:
                branch_lengths = extract_branch_lengths_from_tree(line)
                all_branch_lengths.extend(branch_lengths)

                if branch_lengths:
                    max_bl = max(branch_lengths)
                    mean_bl = np.mean(branch_lengths)

                    tree_stats.append(
                        {
                            "tree_id": tree_idx,
                            "max_branch_length": max_bl,
                            "mean_branch_length": mean_bl,
                            "total_branches": len(branch_lengths),
                        }
                    )

                    # Flag problematic trees (arbitrarily high branch lengths)
                    if max_bl > 1.0:  # Threshold for "high" branch length
                        problematic_trees.append(
                            {
                                "tree_id": tree_idx,
                                "max_branch_length": max_bl,
                                "tree_string": line,
                            }
                        )

    return {
        "all_branch_lengths": all_branch_lengths,
        "tree_stats": tree_stats,
        "problematic_trees": problematic_trees,
    }


def identify_problematic_sequences(tree_file):
    """Identify sequences associated with long branches."""
    if not os.path.exists(tree_file):
        return {}

    sequence_issues = defaultdict(list)

    with open(tree_file, "r") as f:
        for tree_idx, line in enumerate(f, 1):
            line = line.strip()
            if line:
                # Find sequences with very long branches (>1.0)
                long_branch_pattern = r"([A-Z0-9_.]+):([0-9.e-]+)"
                matches = re.findall(long_branch_pattern, line)

                for seq_id, branch_length in matches:
                    bl = float(branch_length)
                    if bl > 1.0:
                        sequence_issues[seq_id].append(
                            {"tree_id": tree_idx, "branch_length": bl}
                        )

    return dict(sequence_issues)


def check_alignment_for_outliers(alignment_file):
    """Check alignment for potential outlier sequences."""
    if not os.path.exists(alignment_file):
        return {}

    sequences = {}
    current_seq = ""
    current_id = ""

    with open(alignment_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_seq:
                    sequences[current_id] = current_seq
                current_id = line[1:].split()[0]  # Get first part of header
                current_seq = ""
            else:
                current_seq += line

        # Add last sequence
        if current_id and current_seq:
            sequences[current_id] = current_seq

    # Analyze sequence properties
    seq_analysis = {}
    seq_lengths = [len(seq) for seq in sequences.values()]
    mean_length = np.mean(seq_lengths)

    for seq_id, seq in sequences.items():
        # Calculate basic stats
        length = len(seq)
        n_count = seq.count("N") + seq.count("n")
        gap_count = seq.count("-")
        gc_content = (
            seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")
        ) / length

        seq_analysis[seq_id] = {
            "length": length,
            "n_count": n_count,
            "gap_count": gap_count,
            "gc_content": gc_content,
            "length_deviation": abs(length - mean_length) / mean_length,
        }

    return seq_analysis


def main():
    """Main analysis function."""
    print("=== HIGH BRANCH LENGTH ANALYSIS ===\n")

    # Files to analyze
    tree_file = "/Users/berksakalli/Projects/automated-window-sliding/results_enhanced_outgroups/best_rooted_trees.newick"
    alignment_file = "/Users/berksakalli/Projects/automated-window-sliding/enhanced_alignment_with_outgroups.fasta"

    # Analyze tree branch lengths
    print("🌳 ANALYZING TREE BRANCH LENGTHS...")
    tree_analysis = analyze_tree_branch_lengths(tree_file)

    if tree_analysis:
        all_bls = tree_analysis["all_branch_lengths"]
        if all_bls:
            print(f"📊 Total branches analyzed: {len(all_bls)}")
            print(f"📏 Branch length statistics:")
            print(f"   Mean: {np.mean(all_bls):.6f}")
            print(f"   Median: {np.median(all_bls):.6f}")
            print(f"   Max: {np.max(all_bls):.6f}")
            print(f"   Min: {np.min(all_bls):.6f}")
            print(f"   95th percentile: {np.percentile(all_bls, 95):.6f}")

            # Count extremely high branch lengths
            very_high = sum(1 for bl in all_bls if bl > 5.0)
            high = sum(1 for bl in all_bls if bl > 1.0)
            print(f"🚨 Branches > 5.0: {very_high}")
            print(f"⚠️  Branches > 1.0: {high}")

        # Show problematic trees
        prob_trees = tree_analysis.get("problematic_trees", [])
        if prob_trees:
            print(f"\n🔴 PROBLEMATIC TREES ({len(prob_trees)} total):")
            for i, tree_info in enumerate(prob_trees[:5]):  # Show first 5
                print(
                    f"   Tree {tree_info['tree_id']}: Max branch = {tree_info['max_branch_length']:.6f}"
                )

            if len(prob_trees) > 5:
                print(f"   ... and {len(prob_trees) - 5} more")

    # Identify problematic sequences
    print(f"\n🧬 IDENTIFYING PROBLEMATIC SEQUENCES...")
    seq_issues = identify_problematic_sequences(tree_file)

    if seq_issues:
        print(f"Found {len(seq_issues)} sequences with long branches:")
        for seq_id, issues in seq_issues.items():
            max_bl = max(issue["branch_length"] for issue in issues)
            tree_count = len(issues)
            print(f"   {seq_id}: Max BL = {max_bl:.6f} (in {tree_count} trees)")

    # Check alignment for outliers
    print(f"\n🔍 ANALYZING ALIGNMENT FOR OUTLIERS...")
    seq_analysis = check_alignment_for_outliers(alignment_file)

    if seq_analysis:
        # Find sequences with unusual properties
        outliers = []
        for seq_id, stats in seq_analysis.items():
            # Flag sequences with high gap content or unusual length
            if (
                stats["gap_count"] / stats["length"] > 0.1  # >10% gaps
                or stats["n_count"] / stats["length"] > 0.05  # >5% Ns
                or stats["length_deviation"] > 0.2
            ):  # >20% length deviation
                outliers.append((seq_id, stats))

        if outliers:
            print(f"Found {len(outliers)} potential alignment outliers:")
            for seq_id, stats in outliers:
                print(f"   {seq_id}:")
                print(
                    f"      Length: {stats['length']} (deviation: {stats['length_deviation']:.3f})"
                )
                print(
                    f"      Gaps: {stats['gap_count']} ({stats['gap_count'] / stats['length'] * 100:.1f}%)"
                )
                print(
                    f"      Ns: {stats['n_count']} ({stats['n_count'] / stats['length'] * 100:.1f}%)"
                )
                print(f"      GC: {stats['gc_content']:.3f}")

    # Recommendations
    print(f"\n=== RECOMMENDATIONS ===")

    # Check if KR074191.1 is the main culprit
    kr074191_issues = seq_issues.get("KR074191.1", [])
    if kr074191_issues:
        max_bl_kr = max(issue["branch_length"] for issue in kr074191_issues)
        print(f"🚨 CRITICAL: KR074191.1 has very long branches (max: {max_bl_kr:.6f})")
        print("   This sequence is likely causing the high branch length issue.")
        print("   Recommendation: Remove KR074191.1 from the alignment")

    # Other recommendations
    if len(seq_issues) > 1:
        print("⚠️  Multiple sequences have long branches:")
        for seq_id in list(seq_issues.keys())[:3]:
            if seq_id != "KR074191.1":
                print(f"   Consider reviewing {seq_id}")

    print("\n🔧 SUGGESTED SOLUTIONS:")
    print("1. Remove highly divergent sequences (especially KR074191.1)")
    print("2. Re-run pipeline with cleaned alignment")
    print("3. Consider using different outgroup strategy")
    print("4. Apply branch length constraints in IQ-TREE")
    print("5. Use molecular clock models if appropriate")

    # Generate cleanup script
    print(f"\n📝 Creating cleanup script...")

    cleanup_script = """#!/usr/bin/env python3
# Script to create cleaned alignment without problematic sequences

from Bio import SeqIO
import sys

def clean_alignment(input_file, output_file, sequences_to_remove):
    '''Remove problematic sequences from alignment.'''
    
    kept_sequences = []
    removed_count = 0
    
    for record in SeqIO.parse(input_file, "fasta"):
        # Check if this sequence should be removed
        should_remove = False
        for seq_to_remove in sequences_to_remove:
            if seq_to_remove in record.id:
                should_remove = True
                removed_count += 1
                print(f"Removing: {record.id}")
                break
        
        if not should_remove:
            kept_sequences.append(record)
    
    # Write cleaned alignment
    SeqIO.write(kept_sequences, output_file, "fasta")
    
    print(f"Kept {len(kept_sequences)} sequences")
    print(f"Removed {removed_count} sequences")
    print(f"Cleaned alignment saved to: {output_file}")

if __name__ == "__main__":
    input_file = "/Users/berksakalli/Projects/automated-window-sliding/enhanced_alignment_with_outgroups.fasta"
    output_file = "/Users/berksakalli/Projects/automated-window-sliding/cleaned_alignment_with_outgroups.fasta"
    
    # Sequences to remove (add/remove as needed)
    sequences_to_remove = ["KR074191.1"]  # Main culprit
    
    clean_alignment(input_file, output_file, sequences_to_remove)
"""

    with open(
        "/Users/berksakalli/Projects/automated-window-sliding/clean_alignment.py", "w"
    ) as f:
        f.write(cleanup_script)

    print("Created: clean_alignment.py")

    print(f"\n=== NEXT STEPS ===")
    print("1. Run: python3 clean_alignment.py")
    print("2. Update parameter file to use cleaned_alignment_with_outgroups.fasta")
    print("3. Re-run pipeline with cleaned alignment")
    print("4. Verify branch lengths are reasonable")


if __name__ == "__main__":
    main()
