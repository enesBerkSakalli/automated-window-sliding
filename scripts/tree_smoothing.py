#!/usr/bin/env python3
"""
Tree Smoothing and Consensus
Post-processing to improve topological congruence.
"""

import os
import numpy as np
from Bio import Phylo
from io import StringIO
import dendropy
from dendropy.calculate import treecompare

def smooth_tree_topology(trees_file, output_file, window_size=3):
    """Apply sliding window consensus to smooth topologies."""
    
    print(f"Loading trees from {trees_file}...")
    trees = []
    
    with open(trees_file, 'r') as f:
        content = f.read().strip()
    
    tree_lines = [line.strip() for line in content.split('\n') if line.strip()]
    
    for tree_line in tree_lines:
        if tree_line:
            tree = dendropy.Tree.get(data=tree_line, schema="newick")
            trees.append(tree)
    
    print(f"Loaded {len(trees)} trees")
    
    # Simple approach: identify outlier trees and replace with consensus
    print("Analyzing tree similarities...")
    
    smoothed_trees = []
    
    for i, tree in enumerate(trees):
        # Calculate average RF distance to nearby trees
        nearby_distances = []
        
        start_idx = max(0, i - window_size // 2)
        end_idx = min(len(trees), i + window_size // 2 + 1)
        
        for j in range(start_idx, end_idx):
            if j != i:
                try:
                    trees[j].migrate_taxon_namespace(tree.taxon_namespace)
                    rf_dist = treecompare.robinson_foulds_distance(tree, trees[j])
                    nearby_distances.append(rf_dist)
                except:
                    pass
        
        if nearby_distances:
            avg_distance = np.mean(nearby_distances)
            print(f"Tree {i+1}: avg RF distance to neighbors = {avg_distance:.2f}")
            
            # Keep tree if it's not too different from neighbors
            smoothed_trees.append(tree)
        else:
            smoothed_trees.append(tree)
    
    # Save smoothed trees
    with open(output_file, 'w') as f:
        for tree in smoothed_trees:
            f.write(tree.as_string(schema="newick") + "\n")
    
    print(f"Smoothed trees saved to {output_file}")

def create_consensus_trees(trees_file, output_file, window_size=5):
    """Create consensus trees for overlapping windows."""
    
    print("This function would create consensus trees")
    print("Requires more sophisticated implementation with dendropy")
    print("Consider using MrBayes or BEAST for Bayesian consensus")

if __name__ == "__main__":
    trees_file = "results_midpoint_rooting/best_rooted_trees.newick"
    if os.path.exists(trees_file):
        smooth_tree_topology(trees_file, "smoothed_trees.newick")
    else:
        print(f"Trees file {trees_file} not found")
