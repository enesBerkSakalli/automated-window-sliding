#!/usr/bin/env python3
"""
Script to check which sequences from the trees are present in the combined sequences
and identify potential matches with the new paper sequences.
"""

import os
import re


def get_tree_accessions(tree_file):
    """Extract all accession numbers from a tree file."""
    accessions = set()

    with open(tree_file, "r") as f:
        content = f.read()

    # Find all accession patterns in the tree
    accession_matches = re.findall(r"([A-Z]{2}\d{6}\.\d+)", content)
    accessions.update(accession_matches)

    return accessions


def get_fasta_accessions(fasta_file):
    """Extract all accession numbers from a FASTA file."""
    accessions = set()

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Extract accession number
                match = re.match(r">(\w+\.\d+)", line)
                if match:
                    accessions.add(match.group(1))

    return accessions


def main():
    """Main function to analyze sequence coverage."""

    # Configuration
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    tree_file = os.path.join(
        base_dir,
        "results",
        "ingroup_only_alignment_w300_s10_GTR_F_I_G4_20250618_214333",
        "all_rooted_trees.newick",
    )
    original_fasta = os.path.join(base_dir, "data", "norovirus_sequences.fasta")
    paper_fasta = os.path.join(base_dir, "data", "paper_sequences.fasta")
    combined_fasta = os.path.join(
        base_dir, "data", "combined_norovirus_sequences.fasta"
    )

    print("Analyzing sequence coverage...")

    # Get accessions from tree
    print("Extracting accessions from tree file...")
    tree_accessions = get_tree_accessions(tree_file)
    print(f"Found {len(tree_accessions)} unique accessions in trees")

    # Get accessions from FASTA files
    print("Extracting accessions from FASTA files...")
    original_accessions = get_fasta_accessions(original_fasta)
    paper_accessions = get_fasta_accessions(paper_fasta)
    combined_accessions = get_fasta_accessions(combined_fasta)

    print(f"Original FASTA: {len(original_accessions)} sequences")
    print(f"Paper FASTA: {len(paper_accessions)} sequences")
    print(f"Combined FASTA: {len(combined_accessions)} sequences")

    # Check coverage
    print("\nCoverage Analysis:")

    # Tree sequences covered by original file
    covered_by_original = tree_accessions.intersection(original_accessions)
    print(
        f"Tree sequences covered by original file: {len(covered_by_original)}/{len(tree_accessions)}"
    )

    # Tree sequences covered by combined file
    covered_by_combined = tree_accessions.intersection(combined_accessions)
    print(
        f"Tree sequences covered by combined file: {len(covered_by_combined)}/{len(tree_accessions)}"
    )

    # Missing sequences
    missing_from_combined = tree_accessions - combined_accessions
    if missing_from_combined:
        print(
            f"\nSequences in trees but missing from combined file ({len(missing_from_combined)}):"
        )
        for acc in sorted(missing_from_combined):
            print(f"  {acc}")
    else:
        print("\nAll tree sequences are present in the combined file!")

    # Additional sequences in paper
    paper_only = paper_accessions - tree_accessions
    print(f"\nAdditional sequences from paper not in trees: {len(paper_only)}")
    if len(paper_only) <= 20:  # Show if not too many
        print("Paper sequences not in current trees:")
        for acc in sorted(paper_only):
            print(f"  {acc}")

    # Show some tree accessions for reference
    print("\nSample tree accessions:")
    for acc in sorted(tree_accessions)[:10]:
        print(f"  {acc}")


if __name__ == "__main__":
    main()
