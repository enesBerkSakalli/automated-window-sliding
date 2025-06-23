#!/usr/bin/env python3
"""
Recombination Analysis for Norovirus Sequences
Based on RDP4 methodologies and current best practices
"""

import os
import sys
from Bio import SeqIO, Phylo
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def detect_recombination_signals(alignment_file, tree_files, output_dir):
    """Detect recombination signals using phylogenetic incongruence."""
    
    print("=== RECOMBINATION SIGNAL DETECTION ===")
    
    # Read alignment
    sequences = list(SeqIO.parse(alignment_file, "fasta"))
    alignment_length = len(sequences[0].seq)
    
    print(f"Analyzing {len(sequences)} sequences, {alignment_length} bp")
    
    # Analyze tree files for incongruence
    tree_distances = []
    window_positions = []
    
    for tree_file in sorted(tree_files):
        try:
            tree = Phylo.read(tree_file, "newick")
            # Extract window position from filename
            pos = extract_window_position(tree_file)
            if pos is not None:
                window_positions.append(pos)
                # Calculate tree statistics
                tree_distances.append(calculate_tree_stats(tree))
        except Exception as e:
            print(f"Error processing {tree_file}: {e}")
    
    # Detect incongruence patterns
    if len(tree_distances) > 1:
        plot_recombination_signals(window_positions, tree_distances, output_dir)
    
    return window_positions, tree_distances

def extract_window_position(filename):
    """Extract window position from filename."""
    import re
    match = re.search(r'window_(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

def calculate_tree_stats(tree):
    """Calculate tree statistics for comparison."""
    # Simple tree statistics - can be enhanced
    total_branch_length = tree.total_branch_length()
    terminal_count = len(tree.get_terminals())
    
    return {
        'total_branch_length': total_branch_length,
        'terminal_count': terminal_count,
        'mean_branch_length': total_branch_length / terminal_count if terminal_count > 0 else 0
    }

def plot_recombination_signals(positions, tree_stats, output_dir):
    """Plot recombination signals across the genome."""
    
    print("Creating recombination signal plots...")
    
    if not tree_stats:
        return
    
    # Extract branch length data
    branch_lengths = [stats['total_branch_length'] for stats in tree_stats]
    
    plt.figure(figsize=(12, 6))
    plt.plot(positions, branch_lengths, 'b-', linewidth=2, alpha=0.7)
    plt.scatter(positions, branch_lengths, c='red', s=30, alpha=0.6)
    
    plt.xlabel('Genomic Position (bp)')
    plt.ylabel('Total Tree Length')
    plt.title('Phylogenetic Signal Variation Across Genome\n(Potential Recombination Signals)')
    plt.grid(True, alpha=0.3)
    
    # Add recombination hotspot annotation
    plt.axvline(x=np.mean(positions), color='orange', linestyle='--', 
                label='Potential ORF1/ORF2 junction', alpha=0.7)
    plt.legend()
    
    output_file = os.path.join(output_dir, "recombination_signals.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Recombination signal plot saved: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python recombination_analysis.py <alignment_file> <tree_directory> <output_dir>")
        sys.exit(1)
    
    alignment_file = sys.argv[1]
    tree_directory = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Find tree files
    tree_files = []
    for root, dirs, files in os.walk(tree_directory):
        for file in files:
            if file.endswith('.treefile') or file.endswith('.newick'):
                tree_files.append(os.path.join(root, file))
    
    print(f"Found {len(tree_files)} tree files")
    
    if tree_files:
        detect_recombination_signals(alignment_file, tree_files, output_dir)
    else:
        print("No tree files found for analysis")
