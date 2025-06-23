#!/usr/bin/env python3
"""
Constrained Phylogenetic Analysis
Uses backbone tree to constrain sliding window analysis.
"""

import subprocess
import os
from Bio import Phylo

def create_backbone_tree():
    """Create backbone tree from full alignment."""
    
    input_file = "ingroup_only_alignment.fasta"
    output_tree = "backbone_tree.newick"
    
    print("Creating backbone tree from full alignment...")
    
    # Use IQ-TREE to create backbone tree
    cmd = [
        "iqtree2",
        "-s", input_file,
        "-m", "GTR+F+I+G4",
        "-pre", "backbone",
        "-redo"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("Backbone tree created successfully")
        
        # Rename output file
        if os.path.exists("backbone.treefile"):
            os.rename("backbone.treefile", output_tree)
            print(f"Backbone tree saved as {output_tree}")
        
        return output_tree
        
    except subprocess.CalledProcessError as e:
        print(f"Error creating backbone tree: {e}")
        return None
    except FileNotFoundError:
        print("IQ-TREE not found. Please install IQ-TREE to create backbone tree.")
        return None

def apply_constraints_to_windows():
    """Apply topological constraints to sliding window analysis."""
    
    print("To apply constraints in IQ-TREE:")
    print("1. Create backbone tree with create_backbone_tree()")
    print("2. Use -g backbone_tree.newick in IQ-TREE commands")
    print("3. This constrains major topology while allowing minor rearrangements")
    
    return "Use -g option in IQ-TREE for constrained analysis"

if __name__ == "__main__":
    backbone_file = create_backbone_tree()
    if backbone_file:
        apply_constraints_to_windows()
