#!/usr/bin/env python3
"""
Manual Constraint-Guided Analysis for Optimized Windows

This script takes the generated sliding windows and applies backbone constraints
followed by RootDigger exhaustive search.
"""

import os
import subprocess
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def run_constraint_analysis():
    """Run constraint-guided analysis on the optimized windows."""

    # Paths
    results_dir = Path("results_optimized_200_20")
    backbone_tree = results_dir / "backbone_analysis" / "backbone_tree.treefile"
    window_alignments = list(
        (results_dir / "sliding_windows" / "alignments").glob("*.fasta")
    )
    constraint_dir = results_dir / "constraint_trees"
    constraint_dir.mkdir(exist_ok=True)

    logger.info(
        f"🔗 Running constraint-guided reconstruction on {len(window_alignments)} windows"
    )
    logger.info(f"🌳 Using backbone tree: {backbone_tree}")

    constraint_results = []

    for i, fasta_file in enumerate(sorted(window_alignments), 1):
        logger.info(
            f"Processing window {i}/{len(window_alignments)}: {fasta_file.name}"
        )

        window_name = fasta_file.stem
        output_prefix = constraint_dir / f"window_{window_name}_constrained"

        # Run IQ-TREE with backbone constraint
        cmd = [
            "iqtree2",
            "-s",
            str(fasta_file),
            "-g",
            str(backbone_tree),  # Backbone constraint
            "-m",
            "MFP",
            "--prefix",
            str(output_prefix),
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            if result.returncode == 0:
                constraint_tree = f"{output_prefix}.treefile"
                if os.path.exists(constraint_tree):
                    constraint_results.append(
                        {
                            "window": window_name,
                            "tree": constraint_tree,
                            "alignment": str(fasta_file),
                            "status": "success",
                        }
                    )
                    logger.info(f"  ✅ Constraint tree generated: {constraint_tree}")
                else:
                    logger.warning(f"  ⚠️ Tree file not found for {fasta_file.name}")
            else:
                logger.warning(
                    f"  ❌ IQ-TREE failed for {fasta_file.name}: {result.stderr}"
                )

        except subprocess.TimeoutExpired:
            logger.warning(f"  ⏰ Timeout for {fasta_file.name}")
        except Exception as e:
            logger.error(f"  ❌ Error processing {fasta_file.name}: {e}")

    logger.info(
        f"✅ Constraint analysis completed: {len(constraint_results)} successful trees"
    )
    return constraint_results


def run_rootdigger_analysis(constraint_results):
    """Apply RootDigger to the constraint trees."""

    results_dir = Path("results_optimized_200_20")
    rooting_dir = results_dir / "rootdigger_exhaustive"
    rooting_dir.mkdir(exist_ok=True)

    logger.info(
        f"🌿 Applying RootDigger exhaustive search to {len(constraint_results)} trees"
    )

    rooted_results = []

    for i, result in enumerate(constraint_results, 1):
        tree_file = result["tree"]
        alignment_file = result["alignment"]
        window_name = result["window"]

        logger.info(f"Processing {i}/{len(constraint_results)}: {window_name}")

        output_prefix = rooting_dir / f"window_{window_name}_rootdigger"

        # Run RootDigger
        cmd = [
            "rootdigger",
            "--msa",
            alignment_file,
            "--tree",
            tree_file,
            "--threads",
            "2",
            "--initial-root-strategy",
            "modified-mad",
            "--exhaustive",
            "--silent",
        ]

        try:
            result_proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=600
            )

            if result_proc.returncode == 0:
                # RootDigger outputs tree(s) to stdout - save the result
                rooted_tree_file = f"{output_prefix}.rooted.tree"
                with open(rooted_tree_file, "w") as f:
                    f.write(result_proc.stdout.strip())

                if (
                    os.path.exists(rooted_tree_file)
                    and os.path.getsize(rooted_tree_file) > 0
                ):
                    rooted_results.append(
                        {
                            "window": window_name,
                            "rooted_tree": rooted_tree_file,
                            "status": "success",
                        }
                    )
                    logger.info(f"  ✅ Rooted tree: {rooted_tree_file}")
                else:
                    logger.warning(
                        f"  ⚠️ Empty or missing rooted tree for {window_name}"
                    )
            else:
                logger.warning(
                    f"  ❌ RootDigger failed for {window_name}: {result_proc.stderr}"
                )

        except subprocess.TimeoutExpired:
            logger.warning(f"  ⏰ RootDigger timeout for {window_name}")
        except Exception as e:
            logger.error(f"  ❌ Error with RootDigger for {window_name}: {e}")

    logger.info(f"✅ RootDigger analysis completed: {len(rooted_results)} rooted trees")
    return rooted_results


def validate_results(rooted_results):
    """Validate the rooted trees."""

    logger.info("🔍 Validating rooted trees...")

    from Bio import Phylo

    valid_trees = 0
    total_trees = len(rooted_results)

    for result in rooted_results:
        tree_file = result["rooted_tree"]
        try:
            tree = Phylo.read(tree_file, "newick")
            root = tree.root
            if root and hasattr(root, "clades") and len(root.clades) == 2:
                valid_trees += 1
            else:
                logger.warning(f"Invalid root structure in {tree_file}")
        except Exception as e:
            logger.warning(f"Error reading {tree_file}: {e}")

    success_rate = valid_trees / total_trees * 100 if total_trees > 0 else 0
    logger.info(
        f"✅ Validation complete: {valid_trees}/{total_trees} ({success_rate:.1f}%) properly rooted"
    )

    return {
        "valid_trees": valid_trees,
        "total_trees": total_trees,
        "success_rate": success_rate,
    }


def main():
    """Run the complete manual analysis."""

    logger.info("🚀 Starting Manual Constraint-Guided Analysis")

    # Step 1: Constraint-guided reconstruction
    constraint_results = run_constraint_analysis()

    if not constraint_results:
        logger.error("❌ No constraint trees generated. Cannot proceed.")
        return

    # Step 2: RootDigger exhaustive search
    rooted_results = run_rootdigger_analysis(constraint_results)

    if not rooted_results:
        logger.error("❌ No trees rooted successfully.")
        return

    # Step 3: Validation
    validation = validate_results(rooted_results)

    # Summary
    logger.info("🎉 Manual Analysis Complete!")
    logger.info(f"📊 Sliding Windows: {len(constraint_results)}")
    logger.info(f"🌳 Constraint Trees: {len(constraint_results)}")
    logger.info(f"🌿 Rooted Trees: {len(rooted_results)}")
    logger.info(f"✅ Success Rate: {validation['success_rate']:.1f}%")

    # Save summary
    results_dir = Path("results_optimized_200_20")
    summary_file = results_dir / "MANUAL_ANALYSIS_SUMMARY.md"

    with open(summary_file, "w") as f:
        f.write(f"""# Manual Constraint-Guided Analysis Summary

## Results Overview

- **Sliding Windows Generated**: {len(constraint_results)}
- **Constraint Trees**: {len(constraint_results)}
- **Rooted Trees**: {len(rooted_results)}
- **Validation Success Rate**: {validation["success_rate"]:.1f}%

## Method

1. **Backbone Constraint Generation**: IQ-TREE with full dataset
2. **Sliding Windows**: 200bp windows, 20bp steps (90% overlap)
3. **Constraint-Guided Reconstruction**: IQ-TREE with backbone constraints
4. **RootDigger Exhaustive Search**: Modified MAD strategy
5. **Validation**: Structural validation of rooted trees

## Optimization Features

- ✅ Optimized window parameters (200bp/20bp)
- ✅ Backbone constraint guidance
- ✅ RootDigger exhaustive search
- ✅ Modified MAD rooting strategy
- ✅ Comprehensive validation

---
*Generated by manual constraint-guided analysis pipeline*
""")

    logger.info(f"📄 Summary saved to: {summary_file}")


if __name__ == "__main__":
    main()
