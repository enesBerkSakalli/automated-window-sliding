#!/usr/bin/env python3
"""
Backbone Midpoint Rooting Script
Applies midpoint rooting to phylogenetic trees with backbone constraint.
"""

import sys
import os
from Bio import Phylo


def midpoint_root_tree(input_tree, output_tree):
    """Apply midpoint rooting to a phylogenetic tree."""
    try:
        tree = Phylo.read(input_tree, "newick")
        tree.root_at_midpoint()
        Phylo.write(tree, output_tree, "newick")
        print(f"Successfully applied midpoint rooting to {input_tree}")
        return True
    except Exception as e:
        print(f"Error during midpoint rooting: {e}")
        return False


def main():
    if len(sys.argv) != 3:
        print("Usage: backbone_midpoint_rooting.py <input_tree> <output_tree>")
        sys.exit(1)

    input_tree = sys.argv[1]
    output_tree = sys.argv[2]

    # Get window name from output file path
    window_name = os.path.splitext(os.path.basename(output_tree))[0].replace(
        "_rooted", ""
    )

    print(f"Window: {window_name}")
    print(f"Input tree: {input_tree}")
    print(f"Output tree: {output_tree}")

    if not os.path.exists(input_tree):
        print(f"Error: Input tree file {input_tree} not found")
        sys.exit(1)

    success = midpoint_root_tree(input_tree, output_tree)

    if success:
        print("Backbone-constrained midpoint rooting completed successfully")
        sys.exit(0)
    else:
        print("ERROR: Midpoint rooting failed, using unrooted tree")
        try:
            with open(input_tree, "r") as infile, open(output_tree, "w") as outfile:
                outfile.write(infile.read())
            print("Fallback copy created")
            sys.exit(1)
        except Exception as e:
            print(f"Error creating fallback: {e}")
            sys.exit(1)


if __name__ == "__main__":
    main()
