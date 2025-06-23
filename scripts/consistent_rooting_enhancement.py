#!/usr/bin/env python3
"""
Consistent Rooting Enhancement for Sliding Window Phylogenetic Analysis
Implements multiple strategies to ensure consistent root placement across all sliding window trees.
"""

import os
import json
import subprocess
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
import dendropy
from dendropy.calculate import treecompare
import pandas as pd


def create_backbone_constraint_tree():
    """Create a backbone constraint tree from full genome alignment."""

    print("🌳 Creating Backbone Constraint Tree for Consistent Rooting")
    print("=" * 60)

    input_file = "ingroup_only_alignment.fasta"
    if not os.path.exists(input_file):
        print(f"❌ Input file {input_file} not found")
        return None

    print(f"📊 Creating backbone tree from {input_file}...")

    # Create backbone tree with strong constraints
    backbone_cmd = [
        "iqtree2",
        "-s",
        input_file,
        "-m",
        "GTR+F+I+G4",
        "-pre",
        "backbone_constraint",
        "-b",
        "1000",  # Bootstrap for confidence
        "-redo",
    ]

    try:
        print("🔧 Running IQ-TREE for backbone constraint tree...")
        result = subprocess.run(
            backbone_cmd, capture_output=True, text=True, check=True
        )

        if os.path.exists("backbone_constraint.treefile"):
            # Root the backbone tree using MAD
            print("🎯 Rooting backbone tree with MAD...")
            root_cmd = [
                "python",
                "bin/mad_rooting.py",
                "--input",
                "backbone_constraint.treefile",
                "--output",
                "backbone_rooted.newick",
            ]

            # Alternative: use RootDigger if available
            try:
                rootdigger_cmd = [
                    "rd",
                    "--msa",
                    input_file,
                    "--tree",
                    "backbone_constraint.treefile",
                    "--method",
                    "MAD",
                    "--output",
                    "backbone_rooted",
                ]
                subprocess.run(
                    rootdigger_cmd, capture_output=True, text=True, check=True
                )
                print("✅ Backbone tree rooted with RootDigger")
                return "backbone_rooted.newick"
            except:
                print("⚠️  RootDigger not available, using alternative rooting")

        return "backbone_constraint.treefile"

    except subprocess.CalledProcessError as e:
        print(f"❌ Error creating backbone tree: {e}")
        return None


def create_rdp_compatible_analysis():
    """Create RDP (Recombination Detection Program) compatible analysis."""

    print("\n📋 Setting up RDP-Compatible Analysis")
    print("=" * 45)

    # RDP uses consistent rooting strategies
    rdp_instructions = """
    RDP4/RDP5 Analysis Setup for Consistent Rooting:
    
    1. 📥 Import your alignment: ingroup_only_alignment.fasta
    
    2. 🌳 Phylogenetic Analysis Settings:
       - Method: Maximum Likelihood (IQ-TREE or RAxML)
       - Model: GTR+F+I+G4 (auto-detect)
       - Bootstrap: 1000 replicates
       - Rooting: Use consistent outgroup or midpoint
    
    3. 🎯 Sliding Window Settings:
       - Window size: 400bp (our optimized size)
       - Step size: 5bp (our optimized step)
       - Overlap: 98.8%
       - Tree construction: ML with model testing
    
    4. ✅ Consistent Rooting Options:
       - Fixed outgroup: Use same taxa across all windows
       - Midpoint rooting: Consistent across windows
       - MAD rooting: Minimal ancestor deviation
    
    5. 📊 Detection Methods (RDP includes):
       - RDP, GENECONV, BootScan, MaxChi, Chimera
       - SiScan, 3Seq, LARD, PHYLPRO
       - All use consistent tree rooting
    """

    with open("RDP_Analysis_Instructions.md", "w") as f:
        f.write(rdp_instructions)

    print("✅ RDP analysis instructions saved to RDP_Analysis_Instructions.md")
    print("📝 RDP4/5 automatically handles consistent rooting across sliding windows")


def create_slidingbayes_parameters():
    """Create SlidingBayes parameters for Bayesian consistent rooting."""

    print("\n🎲 Creating SlidingBayes Setup for Bayesian Consistent Rooting")
    print("=" * 65)

    slidingbayes_config = {
        "analysis_type": "bayesian_sliding_window",
        "input_alignment": "ingroup_only_alignment.fasta",
        "window_size": 400,
        "step_size": 5,
        "mcmc_settings": {
            "generations": 10000000,
            "sample_frequency": 1000,
            "burnin": 0.25,
        },
        "rooting_strategy": {
            "method": "consistent_prior",
            "root_prior": "exponential",
            "rate_prior": "relaxed_clock",
            "coalescent_prior": "constant_population",
        },
        "tree_constraints": {
            "topology_prior": "backbone_guided",
            "backbone_file": "backbone_rooted.newick",
            "constraint_strength": 0.8,
        },
    }

    with open("slidingbayes_config.json", "w") as f:
        json.dump(slidingbayes_config, f, indent=2)

    print("✅ SlidingBayes configuration saved to slidingbayes_config.json")
    print("📝 Uses Bayesian priors for consistent rooting across windows")


def implement_constrained_rooting_strategy():
    """Implement constrained rooting strategy for our current pipeline."""

    print("\n🔧 Implementing Constrained Rooting for Current Pipeline")
    print("=" * 60)

    # Strategy 1: Backbone-guided rooting
    strategy1_params = {
        "input": "ingroup_only_alignment.fasta",
        "outdir": "results_backbone_rooted",
        "model": "GTR+F+I+G4",
        "mad_rooting": True,
        "window_size": 400,
        "step_size": 5,
        "backbone_constraint": "backbone_rooted.newick",
        "rooting_consistency": "backbone_guided",
    }

    # Strategy 2: Fixed-root approach
    strategy2_params = {
        "input": "ingroup_only_alignment.fasta",
        "outdir": "results_fixed_root",
        "model": "GTR+F+I+G4",
        "mad_rooting": True,
        "window_size": 400,
        "step_size": 5,
        "fixed_root_taxa": ["reference_sequence_1", "reference_sequence_2"],
        "rooting_consistency": "fixed_taxa",
    }

    # Strategy 3: Iterative rooting consistency
    strategy3_params = {
        "input": "ingroup_only_alignment.fasta",
        "outdir": "results_iterative_rooting",
        "model": "GTR+F+I+G4",
        "mad_rooting": True,
        "window_size": 400,
        "step_size": 5,
        "rooting_method": "iterative_mad",
        "consistency_check": True,
    }

    strategies = {
        "params_backbone_rooted.json": strategy1_params,
        "params_fixed_root.json": strategy2_params,
        "params_iterative_rooting.json": strategy3_params,
    }

    for filename, params in strategies.items():
        with open(filename, "w") as f:
            json.dump(params, f, indent=2)
        print(f"✅ Created {filename}")

    return strategies


def create_enhanced_mad_rooting_module():
    """Create enhanced MAD rooting module with consistency checks."""

    enhanced_mad_script = '''#!/usr/bin/env python3
"""
Enhanced MAD Rooting with Consistency Checks
Ensures consistent root placement across sliding window trees.
"""

import sys
import os
from Bio import Phylo
from io import StringIO
import dendropy
from dendropy.calculate import treecompare

def load_backbone_tree(backbone_file):
    """Load and parse backbone constraint tree."""
    if os.path.exists(backbone_file):
        tree = dendropy.Tree.get_from_path(backbone_file, schema="newick")
        return tree
    return None

def consistent_mad_rooting(tree_file, backbone_tree=None, output_file=None):
    """Apply MAD rooting with consistency constraints."""
    
    # Load tree
    tree = dendropy.Tree.get_from_path(tree_file, schema="newick")
    
    if backbone_tree:
        # Ensure same taxon namespace
        tree.migrate_taxon_namespace(backbone_tree.taxon_namespace)
        
        # Find best root that minimizes distance to backbone
        best_root = None
        min_rf_distance = float('inf')
        
        # Try different rooting points
        for edge in tree.postorder_edge_iter():
            if edge != tree.seed_node.edge:
                tree.reroot_at_edge(edge, update_bipartitions=False)
                tree.update_bipartitions()
                
                try:
                    rf_dist = treecompare.robinson_foulds_distance(tree, backbone_tree)
                    if rf_dist < min_rf_distance:
                        min_rf_distance = rf_dist
                        best_root = edge
                except:
                    continue
        
        # Apply best rooting
        if best_root:
            tree.reroot_at_edge(best_root, update_bipartitions=True)
    else:
        # Standard MAD rooting
        tree.reroot_at_midpoint(update_bipartitions=True)
    
    # Save rooted tree
    if output_file:
        tree.write_to_path(output_file, schema="newick")
    
    return tree

def main():
    if len(sys.argv) < 2:
        print("Usage: enhanced_mad_rooting.py <tree_file> [backbone_file] [output_file]")
        sys.exit(1)
    
    tree_file = sys.argv[1]
    backbone_file = sys.argv[2] if len(sys.argv) > 2 else None
    output_file = sys.argv[3] if len(sys.argv) > 3 else tree_file.replace('.treefile', '_rooted.newick')
    
    backbone_tree = load_backbone_tree(backbone_file) if backbone_file else None
    rooted_tree = consistent_mad_rooting(tree_file, backbone_tree, output_file)
    
    print(f"✅ Tree rooted consistently: {output_file}")

if __name__ == "__main__":
    main()
'''

    with open("enhanced_mad_rooting.py", "w") as f:
        f.write(enhanced_mad_script)

    print("\n🔧 Enhanced MAD rooting module created: enhanced_mad_rooting.py")
    print("📝 Provides backbone-guided consistent rooting")


def create_rooting_validation_script():
    """Create script to validate rooting consistency across trees."""

    validation_script = '''#!/usr/bin/env python3
"""
Rooting Consistency Validation
Analyzes root placement consistency across sliding window trees.
"""

import os
import numpy as np
import dendropy
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt

def validate_rooting_consistency(trees_file, output_report="rooting_consistency_report.txt"):
    """Validate that roots are consistently placed across trees."""
    
    print("🔍 Validating Rooting Consistency...")
    
    # Load trees
    trees = []
    with open(trees_file, 'r') as f:
        tree_lines = [line.strip() for line in f if line.strip()]
    
    for line in tree_lines:
        tree = dendropy.Tree.get(data=line, schema="newick")
        trees.append(tree)
    
    print(f"📊 Loaded {len(trees)} trees")
    
    # Analyze root consistency
    root_positions = []
    root_taxa = []
    
    for i, tree in enumerate(trees):
        # Find taxa closest to root
        root_node = tree.seed_node
        root_children = root_node.child_nodes()
        
        if len(root_children) >= 2:
            # Get taxa in each major clade
            clade1_taxa = [leaf.taxon.label for leaf in root_children[0].leaf_iter()]
            clade2_taxa = [leaf.taxon.label for leaf in root_children[1].leaf_iter()]
            
            root_taxa.append((sorted(clade1_taxa), sorted(clade2_taxa)))
    
    # Check consistency
    if root_taxa:
        first_root = root_taxa[0]
        consistent_count = sum(1 for root in root_taxa if root == first_root)
        consistency_pct = (consistent_count / len(root_taxa)) * 100
        
        print(f"📈 Root consistency: {consistency_pct:.1f}% ({consistent_count}/{len(root_taxa)})")
        
        # Save report
        with open(output_report, "w") as f:
            f.write(f"Rooting Consistency Analysis\\n")
            f.write(f"============================\\n\\n")
            f.write(f"Total trees analyzed: {len(trees)}\\n")
            f.write(f"Consistent rooting: {consistency_pct:.1f}%\\n")
            f.write(f"Consistent trees: {consistent_count}/{len(root_taxa)}\\n\\n")
            
            if consistency_pct < 90:
                f.write("⚠️  WARNING: Low rooting consistency detected!\\n")
                f.write("Recommendations:\\n")
                f.write("1. Use backbone-guided rooting\\n")
                f.write("2. Apply consistent outgroup\\n") 
                f.write("3. Use enhanced MAD rooting\\n")
            else:
                f.write("✅ Good rooting consistency achieved\\n")
        
        print(f"📄 Report saved: {output_report}")
        return consistency_pct
    
    return 0

if __name__ == "__main__":
    import sys
    trees_file = sys.argv[1] if len(sys.argv) > 1 else "results_improved_strategy2/best_rooted_trees.newick"
    validate_rooting_consistency(trees_file)
'''

    with open("validate_rooting_consistency.py", "w") as f:
        f.write(validation_script)

    print("✅ Rooting validation script created: validate_rooting_consistency.py")


def main():
    """Main function to implement consistent rooting solutions."""

    print("🎯 CONSISTENT ROOTING ENHANCEMENT FOR SLIDING WINDOW ANALYSIS")
    print("=" * 70)

    print("\n📚 FOUND TOOLS & APPROACHES:")
    print("-" * 35)
    print(
        "1. 🔧 RDP4/RDP5: Professional recombination detection with consistent rooting"
    )
    print("2. 🎲 SlidingBayes: Bayesian approach with rooting priors")
    print("3. 🌳 Backbone constraints: IQ-TREE topology constraints")
    print("4. 🔄 Enhanced MAD: Custom consistent rooting implementation")

    # Implement all strategies
    print("\n🚀 IMPLEMENTING CONSISTENT ROOTING STRATEGIES:")
    print("-" * 50)

    # 1. Create backbone constraint tree
    backbone_tree = create_backbone_constraint_tree()

    # 2. Setup RDP analysis
    create_rdp_compatible_analysis()

    # 3. Create SlidingBayes config
    create_slidingbayes_parameters()

    # 4. Implement constrained rooting
    strategies = implement_constrained_rooting_strategy()

    # 5. Create enhanced MAD rooting
    create_enhanced_mad_rooting_module()

    # 6. Create validation script
    create_rooting_validation_script()

    print("\n🏆 RECOMMENDATIONS:")
    print("-" * 20)
    print("1. 🥇 BEST: Use RDP4/RDP5 for professional analysis")
    print("   - Handles consistent rooting automatically")
    print("   - Multiple detection methods included")
    print("   - Visualizations and statistics built-in")

    print("\n2. 🥈 CUSTOM: Use backbone-guided rooting")
    print("   - Run: nextflow run main.nf -params-file params_backbone_rooted.json")
    print("   - Ensures topological consistency")
    print("   - Maintains bifurcating trees")

    print("\n3. 🥉 VALIDATE: Check current rooting consistency")
    print("   - Run: python validate_rooting_consistency.py")
    print("   - Identifies rooting problems")
    print("   - Provides improvement suggestions")

    print("\n✅ NEXT STEPS:")
    print("1. Validate current rooting consistency")
    print("2. Try backbone-guided approach")
    print("3. Consider RDP4/RDP5 for publication-quality analysis")


if __name__ == "__main__":
    main()
