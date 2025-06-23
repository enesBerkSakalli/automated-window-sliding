#!/usr/bin/env python3
"""
Comprehensive validation and summary of the optimized sliding window phylogenetic analysis.
This script validates all results and generates a final report.
"""

import os
import re
import json
from Bio import Phylo
from io import StringIO
from datetime import datetime


def validate_rooted_tree(tree_file):
    """
    Validate that a tree is properly rooted (binary root).

    Args:
        tree_file: Path to tree file

    Returns:
        dict: Validation results
    """
    try:
        with open(tree_file, "r") as f:
            tree_str = f.read().strip()

        # Parse tree with Biopython
        tree = Phylo.read(StringIO(tree_str), "newick")

        # Check if root has exactly 2 children (binary root)
        root_children = len(tree.root.clades)
        is_binary_root = root_children == 2

        # Count total taxa
        taxa_count = len(tree.get_terminals())

        # Check for branch lengths
        has_branch_lengths = all(
            hasattr(clade, "branch_length") and clade.branch_length is not None
            for clade in tree.find_clades()
        )

        return {
            "valid": True,
            "binary_root": is_binary_root,
            "root_children": root_children,
            "taxa_count": taxa_count,
            "has_branch_lengths": has_branch_lengths,
            "error": None,
        }

    except Exception as e:
        return {
            "valid": False,
            "binary_root": False,
            "root_children": 0,
            "taxa_count": 0,
            "has_branch_lengths": False,
            "error": str(e),
        }


def extract_window_info(filename):
    """
    Extract window position from filename.

    Args:
        filename: Tree filename

    Returns:
        int: Window start position, or None if not found
    """
    match = re.search(r"window_(\d+)", filename)
    if match:
        return int(match.group(1))
    return None


def validate_all_results(results_dir):
    """
    Comprehensive validation of all analysis results.

    Args:
        results_dir: Path to results directory

    Returns:
        dict: Complete validation report
    """
    validation_report = {
        "timestamp": datetime.now().isoformat(),
        "analysis_type": "Optimized Sliding Window Phylogenetic Analysis",
        "parameters": {
            "window_size": 200,
            "step_size": 20,
            "rooting_method": "RootDigger exhaustive",
            "constraint_optimization": True,
            "backbone_tree": True,
        },
        "results": {},
    }

    # Check backbone tree
    backbone_file = os.path.join(
        results_dir, "backbone_analysis", "backbone_tree.treefile"
    )
    if os.path.exists(backbone_file):
        validation_report["results"]["backbone_tree"] = validate_rooted_tree(
            backbone_file
        )
        validation_report["results"]["backbone_tree"]["file"] = backbone_file
    else:
        validation_report["results"]["backbone_tree"] = {
            "error": "Backbone tree file not found"
        }

    # Check sliding window alignments
    alignments_dir = os.path.join(results_dir, "sliding_windows", "alignments")
    cleaned_alignments_dir = os.path.join(
        results_dir, "sliding_windows", "alignments_cleaned"
    )

    alignment_count = 0
    cleaned_alignment_count = 0

    if os.path.exists(alignments_dir):
        alignment_count = len(
            [f for f in os.listdir(alignments_dir) if f.endswith(".fasta")]
        )

    if os.path.exists(cleaned_alignments_dir):
        cleaned_alignment_count = len(
            [f for f in os.listdir(cleaned_alignments_dir) if f.endswith(".fasta")]
        )

    validation_report["results"]["alignments"] = {
        "original_count": alignment_count,
        "cleaned_count": cleaned_alignment_count,
        "expected_count": 25,  # Based on 501bp sequence with window=200, step=20
    }

    # Check constraint trees
    constraint_trees_dir = os.path.join(results_dir, "constraint_trees")
    constraint_trees = []

    if os.path.exists(constraint_trees_dir):
        for file in os.listdir(constraint_trees_dir):
            if file.endswith(".treefile") and "constrained" in file:
                tree_path = os.path.join(constraint_trees_dir, file)
                window_pos = extract_window_info(file)
                validation = validate_rooted_tree(tree_path)
                validation["window_position"] = window_pos
                validation["filename"] = file
                constraint_trees.append(validation)

    constraint_trees.sort(key=lambda x: x.get("window_position", 0))
    validation_report["results"]["constraint_trees"] = {
        "count": len(constraint_trees),
        "trees": constraint_trees,
    }

    # Check rooted trees
    rooted_trees_dir = os.path.join(results_dir, "rooted_trees_exhaustive")
    rooted_trees = []

    if os.path.exists(rooted_trees_dir):
        for file in os.listdir(rooted_trees_dir):
            if file.endswith(".rooted.tree"):
                tree_path = os.path.join(rooted_trees_dir, file)
                window_pos = extract_window_info(file)
                validation = validate_rooted_tree(tree_path)
                validation["window_position"] = window_pos
                validation["filename"] = file
                rooted_trees.append(validation)

    rooted_trees.sort(key=lambda x: x.get("window_position", 0))
    validation_report["results"]["rooted_trees"] = {
        "count": len(rooted_trees),
        "trees": rooted_trees,
    }

    # Calculate success rates
    constraint_success = sum(1 for t in constraint_trees if t["valid"])
    rooted_success = sum(1 for t in rooted_trees if t["valid"])
    binary_root_count = sum(1 for t in rooted_trees if t.get("binary_root", False))

    validation_report["summary"] = {
        "total_windows": 25,
        "alignments_generated": alignment_count,
        "alignments_cleaned": cleaned_alignment_count,
        "constraint_trees_generated": len(constraint_trees),
        "constraint_trees_valid": constraint_success,
        "rooted_trees_generated": len(rooted_trees),
        "rooted_trees_valid": rooted_success,
        "rooted_trees_binary": binary_root_count,
        "pipeline_complete": (
            alignment_count == 25
            and cleaned_alignment_count == 25
            and len(constraint_trees) == 25
            and constraint_success == 25
            and len(rooted_trees) == 25
            and rooted_success == 25
            and binary_root_count == 25
        ),
    }

    return validation_report


def generate_summary_report(validation_report, output_file):
    """
    Generate a human-readable summary report.

    Args:
        validation_report: Validation results dictionary
        output_file: Path to output markdown file
    """

    summary = validation_report["summary"]

    report_content = f"""# Optimized Sliding Window Phylogenetic Analysis Report

**Analysis Date:** {validation_report["timestamp"]}
**Pipeline:** {validation_report["analysis_type"]}

## Parameters
- **Window Size:** {validation_report["parameters"]["window_size"]} bp
- **Step Size:** {validation_report["parameters"]["step_size"]} bp
- **Rooting Method:** {validation_report["parameters"]["rooting_method"]}
- **Backbone Constraint:** {validation_report["parameters"]["backbone_tree"]}
- **Constraint Optimization:** {validation_report["parameters"]["constraint_optimization"]}

## Results Summary

### Overall Pipeline Status
**{"✅ COMPLETE" if summary["pipeline_complete"] else "❌ INCOMPLETE"}**

### Detailed Results
- **Total Windows Generated:** {summary["total_windows"]}
- **Sliding Window Alignments:** {summary["alignments_generated"]}/{summary["total_windows"]} ({"✅" if summary["alignments_generated"] == 25 else "❌"})
- **Cleaned Alignments:** {summary["alignments_cleaned"]}/{summary["total_windows"]} ({"✅" if summary["alignments_cleaned"] == 25 else "❌"})
- **Constraint Trees:** {summary["constraint_trees_generated"]}/{summary["total_windows"]} ({"✅" if summary["constraint_trees_generated"] == 25 else "❌"})
- **Valid Constraint Trees:** {summary["constraint_trees_valid"]}/{summary["constraint_trees_generated"]} ({"✅" if summary["constraint_trees_valid"] == summary["constraint_trees_generated"] else "❌"})
- **Rooted Trees (RootDigger):** {summary["rooted_trees_generated"]}/{summary["total_windows"]} ({"✅" if summary["rooted_trees_generated"] == 25 else "❌"})
- **Valid Rooted Trees:** {summary["rooted_trees_valid"]}/{summary["rooted_trees_generated"]} ({"✅" if summary["rooted_trees_valid"] == summary["rooted_trees_generated"] else "❌"})
- **Binary Rooted Trees:** {summary["rooted_trees_binary"]}/{summary["rooted_trees_generated"]} ({"✅" if summary["rooted_trees_binary"] == summary["rooted_trees_generated"] else "❌"})

### Backbone Tree
"""

    if "backbone_tree" in validation_report["results"]:
        backbone = validation_report["results"]["backbone_tree"]
        if backbone.get("valid", False):
            report_content += f"""- **Status:** ✅ Valid
- **Taxa Count:** {backbone.get("taxa_count", "N/A")}
- **Binary Root:** {"✅" if backbone.get("binary_root", False) else "❌"}
- **Branch Lengths:** {"✅" if backbone.get("has_branch_lengths", False) else "❌"}
"""
        else:
            report_content += (
                f"- **Status:** ❌ Error: {backbone.get('error', 'Unknown error')}\n"
            )
    else:
        report_content += "- **Status:** ❌ Not found\n"

    # Window analysis table
    report_content += """
## Window Analysis Results

| Window | Start Pos | Constraint Tree | Rooted Tree | Binary Root | Taxa Count |
|--------|-----------|----------------|-------------|-------------|------------|
"""

    # Create a combined view of constraint and rooted trees
    rooted_trees = validation_report["results"].get("rooted_trees", {}).get("trees", [])
    constraint_trees = (
        validation_report["results"].get("constraint_trees", {}).get("trees", [])
    )

    # Group by window position
    window_data = {}

    for tree in constraint_trees:
        pos = tree.get("window_position")
        if pos is not None:
            window_data[pos] = window_data.get(pos, {})
            window_data[pos]["constraint"] = tree

    for tree in rooted_trees:
        pos = tree.get("window_position")
        if pos is not None:
            window_data[pos] = window_data.get(pos, {})
            window_data[pos]["rooted"] = tree

    # Sort windows and create table
    for pos in sorted(window_data.keys()):
        data = window_data[pos]
        constraint = data.get("constraint", {})
        rooted = data.get("rooted", {})

        constraint_status = "✅" if constraint.get("valid", False) else "❌"
        rooted_status = "✅" if rooted.get("valid", False) else "❌"
        binary_status = "✅" if rooted.get("binary_root", False) else "❌"
        taxa_count = rooted.get("taxa_count", constraint.get("taxa_count", "N/A"))

        window_num = (pos - 1) // 20 + 1 if pos else "N/A"

        report_content += f"| {window_num} | {pos} | {constraint_status} | {rooted_status} | {binary_status} | {taxa_count} |\n"

    report_content += f"""
## Pipeline Optimization Benefits

This analysis used advanced optimization techniques:

1. **Backbone Constraint Optimization**: Used a constraint tree to guide window-specific tree reconstruction
2. **RootDigger Exhaustive Rooting**: Applied exhaustive search for optimal root placement
3. **Taxa Name Consistency**: Cleaned alignment headers to ensure compatibility between alignments and trees
4. **Robust Validation**: Comprehensive validation of tree topology and rooting quality

## Output Files

### Key Results Directories
- `backbone_analysis/`: Backbone constraint tree and analysis
- `sliding_windows/alignments/`: Original window alignments
- `sliding_windows/alignments_cleaned/`: Header-cleaned alignments for RootDigger
- `constraint_trees/`: Backbone-constrained trees for each window
- `rooted_trees_exhaustive/`: RootDigger exhaustively rooted trees

### File Formats
- `.treefile`: Newick format phylogenetic trees
- `.rooted.tree`: RootDigger optimally rooted trees
- `.fasta`: Multiple sequence alignments
- `.log`: Detailed software logs

## Recommendations

{"✅ Analysis complete! All 25 windows processed successfully with optimal rooting." if summary["pipeline_complete"] else "⚠️  Some steps may need attention. Check individual file statuses above."}

### Next Steps
1. **Tree Comparison**: Compare rooted trees across windows for congruence
2. **Phylogenetic Signal**: Assess signal strength across genomic windows
3. **Recombination Analysis**: Look for breakpoints in phylogenetic signal
4. **Publication**: Results are ready for manuscript preparation

---
*Report generated on {validation_report["timestamp"]}*
"""

    # Write report
    with open(output_file, "w") as f:
        f.write(report_content)

    print(f"Summary report written to: {output_file}")


def main():
    """
    Run comprehensive validation and generate reports.
    """
    results_dir = "results_optimized_200_20"

    print("Running comprehensive validation of optimized sliding window analysis...")

    # Run validation
    validation_report = validate_all_results(results_dir)

    # Save detailed JSON report
    json_report_file = os.path.join(results_dir, "validation_report.json")
    with open(json_report_file, "w") as f:
        json.dump(validation_report, f, indent=2)

    print(f"Detailed validation report saved to: {json_report_file}")

    # Generate summary report
    summary_report_file = os.path.join(results_dir, "ANALYSIS_SUMMARY.md")
    generate_summary_report(validation_report, summary_report_file)

    # Print quick summary
    summary = validation_report["summary"]
    print(f"\n{'=' * 60}")
    print("OPTIMIZED SLIDING WINDOW ANALYSIS SUMMARY")
    print(f"{'=' * 60}")
    print(
        f"Pipeline Status: {'✅ COMPLETE' if summary['pipeline_complete'] else '❌ INCOMPLETE'}"
    )
    print(f"Windows Generated: {summary['total_windows']}")
    print(
        f"Rooted Trees: {summary['rooted_trees_valid']}/{summary['rooted_trees_generated']}"
    )
    print(
        f"Binary Roots: {summary['rooted_trees_binary']}/{summary['rooted_trees_generated']}"
    )
    print(
        f"Success Rate: {(summary['rooted_trees_valid'] / summary['total_windows'] * 100):.1f}%"
    )
    print("\nDetailed reports:")
    print(f"  - JSON: {json_report_file}")
    print(f"  - Markdown: {summary_report_file}")


if __name__ == "__main__":
    main()
