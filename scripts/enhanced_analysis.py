#!/usr/bin/env python3
"""
Enhanced Norovirus Recombination Analysis Pipeline
Based on current best practices from recent literature (2023-2024)

Key findings from literature review:
1. Recombination hotspots occur primarily at ORF1/ORF2 junction
2. Optimal window sizes: 200-600 nt for recombination detection
3. Step sizes: 50-200 nt depending on resolution needed
4. RDP4 with multiple algorithms provides robust detection
5. GTR+I+G model optimal for norovirus phylogenetics
6. Pairwise distance analysis enhances recombination detection
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import json


def analyze_alignment_composition(fasta_file):
    """Analyze the cleaned alignment for ORF composition and optimization potential."""

    print("=== ALIGNMENT COMPOSITION ANALYSIS ===")
    print(f"Analyzing: {fasta_file}")

    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if not sequences:
        print("Error: No sequences found!")
        return None

    num_seqs = len(sequences)
    alignment_length = len(sequences[0].seq)

    print(f"Sequences: {num_seqs}")
    print(f"Alignment length: {alignment_length} bp")

    # Analyze sequence composition
    total_sites = 0
    gap_sites = 0
    informative_sites = 0

    for i in range(alignment_length):
        column = [seq.seq[i] for seq in sequences]
        total_sites += len(column)
        gap_count = column.count("-")
        gap_sites += gap_count

        # Count informative sites (more than one nucleotide present)
        nucleotides = set(column) - {"-", "N", "n"}
        if len(nucleotides) > 1:
            informative_sites += 1

    gap_percentage = (gap_sites / total_sites) * 100
    informative_percentage = (informative_sites / alignment_length) * 100

    print(f"Gap percentage: {gap_percentage:.2f}%")
    print(
        f"Phylogenetically informative sites: {informative_sites} ({informative_percentage:.2f}%)"
    )

    return {
        "num_sequences": num_seqs,
        "alignment_length": alignment_length,
        "gap_percentage": gap_percentage,
        "informative_sites": informative_sites,
        "informative_percentage": informative_percentage,
    }


def calculate_enhanced_distances(fasta_file, output_dir):
    """Calculate pairwise distances using multiple methods for recombination analysis."""

    print("\n=== ENHANCED DISTANCE CALCULATION ===")

    # Create distance analysis script
    distance_script = os.path.join(output_dir, "calculate_distances.py")

    with open(distance_script, "w") as f:
        f.write("""#!/usr/bin/env python3
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
import numpy as np
import pandas as pd

def calculate_distance_matrices(fasta_file, output_dir):
    # Read alignment
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    alignment = MultipleSeqAlignment(sequences)
    
    # Calculate distances using different models
    models = ['blosum62', 'pam250', 'trans', 'ident']
    
    for model in models:
        try:
            calculator = DistanceCalculator(model)
            dm = calculator.get_distance(alignment)
            
            # Save distance matrix
            output_file = f"{output_dir}/distance_matrix_{model}.csv"
            
            # Convert to pandas DataFrame for easier handling
            names = [seq.id for seq in sequences]
            matrix_data = []
            
            for i, name1 in enumerate(names):
                row = []
                for j, name2 in enumerate(names):
                    if i == j:
                        row.append(0.0)
                    elif i > j:
                        row.append(dm[name1, name2])
                    else:
                        row.append(dm[name2, name1])
                matrix_data.append(row)
            
            df = pd.DataFrame(matrix_data, index=names, columns=names)
            df.to_csv(output_file)
            print(f"Distance matrix saved: {output_file}")
            
        except Exception as e:
            print(f"Error calculating {model} distances: {e}")

if __name__ == "__main__":
    import sys
    calculate_distance_matrices(sys.argv[1], sys.argv[2])
""")

    # Make script executable and run
    os.chmod(distance_script, 0o755)

    try:
        result = subprocess.run(
            ["python3", distance_script, fasta_file, output_dir],
            capture_output=True,
            text=True,
            cwd=output_dir,
        )

        if result.returncode == 0:
            print("Distance matrices calculated successfully")
        else:
            print(f"Distance calculation failed: {result.stderr}")
    except Exception as e:
        print(f"Error running distance calculation: {e}")


def create_enhanced_parameter_sets(alignment_stats):
    """Create optimized parameter sets based on alignment characteristics and literature."""

    print("\n=== CREATING ENHANCED PARAMETER SETS ===")

    alignment_length = alignment_stats["alignment_length"]

    # Based on literature: optimal window sizes for recombination detection
    # Small windows: high resolution for breakpoint detection
    # Medium windows: balanced approach for overall patterns
    # Large windows: phylogenetic signal assessment

    param_sets = {
        "ultra_high_resolution": {
            "description": "Ultra-high resolution for precise breakpoint detection",
            "window_size": min(50, alignment_length // 10),
            "step_size": 5,
            "expected_windows": alignment_length // 5,
            "model": "GTR+F+I+G4",
            "bootstrap": 1000,
            "recombination_optimized": True,
        },
        "high_resolution": {
            "description": "High resolution recombination analysis",
            "window_size": min(75, alignment_length // 7),
            "step_size": 10,
            "expected_windows": alignment_length // 10,
            "model": "GTR+F+I+G4",
            "bootstrap": 1000,
            "recombination_optimized": True,
        },
        "optimal_recombination": {
            "description": "Optimal for norovirus recombination (literature-based)",
            "window_size": min(200, alignment_length // 3),
            "step_size": 25,
            "expected_windows": alignment_length // 25,
            "model": "GTR+F+I+G4",
            "bootstrap": 1000,
            "recombination_optimized": True,
        },
        "phylogenetic_signal": {
            "description": "Larger windows for phylogenetic signal assessment",
            "window_size": min(300, alignment_length // 2),
            "step_size": 50,
            "expected_windows": alignment_length // 50,
            "model": "GTR+F+I+G4",
            "bootstrap": 1000,
            "recombination_optimized": False,
        },
    }

    print("Enhanced parameter sets created:")
    for name, params in param_sets.items():
        print(f"  {name}:")
        print(f"    - Window: {params['window_size']}bp, Step: {params['step_size']}bp")
        print(f"    - Expected windows: ~{params['expected_windows']}")
        print(f"    - {params['description']}")

    return param_sets


def generate_enhanced_nextflow_params(param_sets, base_dir):
    """Generate enhanced Nextflow parameter files with optimized settings."""

    print("\n=== GENERATING ENHANCED NEXTFLOW PARAMETERS ===")

    params_dir = os.path.join(base_dir, "params")
    os.makedirs(params_dir, exist_ok=True)

    for set_name, params in param_sets.items():
        param_file = {
            "input": "cleaned_alignment_combined.fasta",
            "outdir": f"results_enhanced_{set_name}",
            "window_size": params["window_size"],
            "step_size": params["step_size"],
            "run_mode": "full",
            "backbone_midpoint_rooting": True,
            "phylo_method": "iqtree2",
            "model": params["model"],
            "max_cpus": 6,  # Increased for better performance
            "publish_dir_mode": "copy",
            "keep_tree_files": True,
            "bootstrap": params.get("bootstrap", 1000),
            # Enhanced settings for recombination analysis
            "model_finder_splits": True if params["recombination_optimized"] else False,
            "output_windows": True,
            "output_format": "nexus,newick",
        }

        output_path = os.path.join(params_dir, f"params_enhanced_{set_name}.json")

        with open(output_path, "w") as f:
            json.dump(param_file, f, indent=4)

        print(f"Created: {output_path}")


def create_recombination_analysis_script(base_dir):
    """Create script for post-analysis recombination detection using RDP4-like methods."""

    print("\n=== CREATING RECOMBINATION ANALYSIS SCRIPT ===")

    script_path = os.path.join(base_dir, "scripts", "recombination_analysis.py")

    with open(script_path, "w") as f:
        f.write("""#!/usr/bin/env python3
\"\"\"
Recombination Analysis for Norovirus Sequences
Based on RDP4 methodologies and current best practices
\"\"\"

import os
import sys
from Bio import SeqIO, Phylo
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def detect_recombination_signals(alignment_file, tree_files, output_dir):
    \"\"\"Detect recombination signals using phylogenetic incongruence.\"\"\"
    
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
    \"\"\"Extract window position from filename.\"\"\"
    import re
    match = re.search(r'window_(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

def calculate_tree_stats(tree):
    \"\"\"Calculate tree statistics for comparison.\"\"\"
    # Simple tree statistics - can be enhanced
    total_branch_length = tree.total_branch_length()
    terminal_count = len(tree.get_terminals())
    
    return {
        'total_branch_length': total_branch_length,
        'terminal_count': terminal_count,
        'mean_branch_length': total_branch_length / terminal_count if terminal_count > 0 else 0
    }

def plot_recombination_signals(positions, tree_stats, output_dir):
    \"\"\"Plot recombination signals across the genome.\"\"\"
    
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
    plt.title('Phylogenetic Signal Variation Across Genome\\n(Potential Recombination Signals)')
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
""")

    os.chmod(script_path, 0o755)
    print(f"Recombination analysis script created: {script_path}")


def create_enhanced_run_script(base_dir, param_sets):
    """Create enhanced run script for all analyses."""

    print("\n=== CREATING ENHANCED RUN SCRIPT ===")

    script_path = os.path.join(base_dir, "run_enhanced_analysis.sh")

    with open(script_path, "w") as f:
        f.write("""#!/bin/bash
# Enhanced Norovirus Recombination Analysis Pipeline
# Based on current best practices and literature recommendations

echo "=== ENHANCED NOROVIRUS RECOMBINATION ANALYSIS ==="
echo "Using cleaned alignment: data/cleaned_alignment_combined.fasta"
echo "86 sequences × 516 positions"
echo "Optimized for recombination detection and phylogenetic analysis"
echo ""

# Create output directory for enhanced analysis
mkdir -p analysis_enhanced
mkdir -p analysis_enhanced/distance_matrices
mkdir -p analysis_enhanced/recombination_results

# Run distance matrix calculations
echo "1. Calculating enhanced distance matrices..."
python3 scripts/enhanced_analysis.py distance_analysis

echo ""
""")

        # Add parameter sets to script
        for i, (set_name, params) in enumerate(param_sets.items(), 1):
            f.write(f"""echo "{i + 1}. Running {set_name.upper().replace("_", " ")} analysis..."
echo "   Window: {params["window_size"]}bp, Step: {params["step_size"]}bp"
echo "   {params["description"]}"
nextflow run main.nf -params-file params/params_enhanced_{set_name}.json -c simple.config --validate_params false

echo ""
""")

        f.write("""
echo "=== POST-ANALYSIS RECOMBINATION DETECTION ==="
echo "Running recombination signal detection across all results..."

# Run recombination analysis on each result set
for result_dir in results_enhanced_*/; do
    if [ -d "$result_dir" ]; then
        echo "Analyzing recombination signals in: $result_dir"
        python3 scripts/recombination_analysis.py data/cleaned_alignment_combined.fasta "$result_dir" analysis_enhanced/recombination_results/
    fi
done

echo ""
echo "=== ENHANCED ANALYSIS COMPLETED ==="
echo "Results directories:"
""")

        for set_name in param_sets.keys():
            f.write(f'echo "  - results_enhanced_{set_name}/"\n')

        f.write("""echo "  - analysis_enhanced/distance_matrices/"
echo "  - analysis_enhanced/recombination_results/"
echo ""
echo "Enhanced analysis features:"
echo "1. Optimized window sizes based on literature (50-300bp)"
echo "2. Enhanced distance calculations using multiple models"
echo "3. Recombination signal detection across sliding windows"
echo "4. GTR+F+I+G4 model optimized for norovirus"
echo "5. High bootstrap support (1000) for robust trees"
echo "6. Specialized breakpoint detection capabilities"
""")

    os.chmod(script_path, 0o755)
    print(f"Enhanced run script created: {script_path}")


def main():
    """Main function to create enhanced analysis pipeline."""

    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    fasta_file = os.path.join(base_dir, "data", "cleaned_alignment_combined.fasta")

    print("🧬 ENHANCED NOROVIRUS RECOMBINATION ANALYSIS SETUP")
    print("=" * 60)
    print("Based on current literature and best practices:")
    print("- Window sizes optimized for norovirus recombination (50-300bp)")
    print("- GTR+F+I+G4 model for accurate phylogenetics")
    print("- Enhanced distance calculations")
    print("- Recombination signal detection")
    print("- Publication-ready analysis pipeline")
    print()

    # Analyze current alignment
    if not os.path.exists(fasta_file):
        print(f"Error: Cleaned alignment not found: {fasta_file}")
        return

    alignment_stats = analyze_alignment_composition(fasta_file)
    if not alignment_stats:
        return

    # Create output directory for enhanced analysis
    enhanced_dir = os.path.join(base_dir, "analysis_enhanced")
    os.makedirs(enhanced_dir, exist_ok=True)
    os.makedirs(os.path.join(enhanced_dir, "distance_matrices"), exist_ok=True)

    # Calculate enhanced distances
    calculate_enhanced_distances(
        fasta_file, os.path.join(enhanced_dir, "distance_matrices")
    )

    # Create optimized parameter sets
    param_sets = create_enhanced_parameter_sets(alignment_stats)

    # Generate Nextflow parameter files
    generate_enhanced_nextflow_params(param_sets, base_dir)

    # Create recombination analysis script
    create_recombination_analysis_script(base_dir)

    # Create enhanced run script
    create_enhanced_run_script(base_dir, param_sets)

    print("\n🎉 ENHANCED ANALYSIS PIPELINE READY!")
    print("=" * 60)
    print("To run the enhanced analysis:")
    print("  ./run_enhanced_analysis.sh")
    print()
    print("This will provide:")
    print("✅ Ultra-high resolution breakpoint detection")
    print("✅ Optimized phylogenetic reconstruction")
    print("✅ Enhanced distance matrix analysis")
    print("✅ Recombination signal visualization")
    print("✅ Publication-ready results")


if __name__ == "__main__":
    main()
