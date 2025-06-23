#!/usr/bin/env python3
"""
Immediate Implementation: IQ-TREE Constraint-Based Rooting for Sliding Window Analysis

This script implements the most promising strategy discovered in the 2024-2025 survey:
using IQ-TREE constraint trees to ensure consistent rooting across sliding windows.

Usage:
    python implement_constraint_rooting.py --alignment ingroup_only_alignment.fasta

Author: Enhanced Norovirus Pipeline
Date: January 2025
"""

import os
import sys
import subprocess
import json
import logging
from pathlib import Path
from typing import List, Dict
import pandas as pd
from Bio import SeqIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("constraint_rooting_implementation.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


class ConstraintRootingImplementation:
    """
    Immediate implementation of constraint-based rooting using IQ-TREE.
    """

    def __init__(self, results_dir: str = "results_constraint_rooting"):
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(exist_ok=True)

        # Check IQ-TREE availability
        if not self._check_iqtree():
            raise RuntimeError(
                "IQ-TREE2 is required but not found. Please install IQ-TREE2."
            )

        logger.info(
            f"Constraint rooting implementation initialized. Results: {self.results_dir}"
        )

    def _check_iqtree(self) -> bool:
        """Check if IQ-TREE2 is available."""
        try:
            result = subprocess.run(
                ["iqtree2", "--version"], capture_output=True, text=True
            )
            return result.returncode == 0
        except FileNotFoundError:
            return False

    def create_backbone_constraint_tree(self, alignment_file: str) -> str:
        """
        Create a high-quality backbone constraint tree from the full alignment.

        Args:
            alignment_file: Path to the full alignment file

        Returns:
            Path to the constraint tree file
        """
        logger.info("Creating backbone constraint tree...")

        backbone_prefix = self.results_dir / "backbone_constraint"

        # Use IQ-TREE2 with model selection and bootstrap for robust backbone
        cmd = [
            "iqtree2",
            "-s",
            alignment_file,
            "-m",
            "MFP",  # Model Finder Plus for best model selection
            "-B",
            "1000",  # 1000 bootstrap replicates
            "-T",
            "AUTO",  # Automatic thread detection
            "--prefix",
            str(backbone_prefix),
            "--quiet",  # Reduce output verbosity
            "-redo",  # Overwrite any existing files
        ]

        logger.info(f"Running IQ-TREE2 for backbone tree: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=3600
            )  # 1 hour timeout

            if result.returncode != 0:
                logger.error(f"IQ-TREE2 failed: {result.stderr}")
                raise RuntimeError("Backbone tree creation failed")

            constraint_tree = f"{backbone_prefix}.treefile"

            if not os.path.exists(constraint_tree):
                raise RuntimeError(
                    f"Expected constraint tree not found: {constraint_tree}"
                )

            logger.info(f"Backbone constraint tree created: {constraint_tree}")
            return constraint_tree

        except subprocess.TimeoutExpired:
            logger.error("IQ-TREE2 timed out after 1 hour")
            raise RuntimeError("Backbone tree creation timed out")

    def generate_sliding_windows(
        self, alignment_file: str, window_size: int = 400, step_size: int = 5
    ) -> List[str]:
        """
        Generate sliding window alignments with optimized parameters.

        Args:
            alignment_file: Path to the alignment file
            window_size: Size of each window (default: 400bp based on previous optimization)
            step_size: Step size between windows (default: 5bp for high overlap)

        Returns:
            List of paths to window alignment files
        """
        logger.info(
            f"Generating sliding windows (window={window_size}, step={step_size})..."
        )

        # Read alignment
        alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, "fasta"))
        seq_length = len(list(alignment.values())[0])

        logger.info(f"Alignment length: {seq_length} bp")

        window_files = []
        window_dir = self.results_dir / "sliding_windows"
        window_dir.mkdir(exist_ok=True)

        window_num = 1
        for start in range(0, seq_length - window_size + 1, step_size):
            end = start + window_size

            # Extract window sequences
            window_alignment = []
            valid_sequences = 0

            for seq_id, seq_record in alignment.items():
                window_seq = str(seq_record.seq)[start:end]

                # Quality check: skip windows with too many gaps or ambiguous bases
                gap_ratio = window_seq.count("-") / len(window_seq)
                n_ratio = window_seq.upper().count("N") / len(window_seq)

                if (
                    gap_ratio < 0.5 and n_ratio < 0.1
                ):  # Less than 50% gaps, less than 10% Ns
                    window_alignment.append(f">{seq_id}\n{window_seq}")
                    valid_sequences += 1

            # Only create window if we have enough sequences
            if valid_sequences >= 4:  # Minimum for meaningful phylogeny
                window_file = (
                    window_dir / f"window_{window_num:03d}_pos_{start}_{end}.fasta"
                )
                with open(window_file, "w") as f:
                    f.write("\n".join(window_alignment))

                window_files.append(str(window_file))
                window_num += 1

        logger.info(f"Generated {len(window_files)} high-quality sliding windows")
        return window_files

    def infer_constraint_guided_trees(
        self, window_files: List[str], constraint_tree: str
    ) -> List[str]:
        """
        Infer phylogenetic trees for each window using the backbone constraint tree.

        Args:
            window_files: List of window alignment files
            constraint_tree: Path to the constraint tree

        Returns:
            List of paths to inferred tree files
        """
        logger.info(
            f"Inferring constraint-guided trees for {len(window_files)} windows..."
        )

        tree_files = []
        failed_windows = []

        for i, window_file in enumerate(window_files, 1):
            window_name = Path(window_file).stem
            output_prefix = self.results_dir / f"{window_name}_constrained"

            cmd = [
                "iqtree2",
                "-s",
                window_file,
                "-g",
                constraint_tree,  # Critical: use constraint tree
                "-m",
                "GTR+I+G",  # Use a reasonable model for speed
                "--prefix",
                str(output_prefix),
                "--quiet",
            ]

            try:
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=600
                )  # 10-minute timeout per window

                if result.returncode == 0:
                    tree_file = f"{output_prefix}.treefile"
                    if os.path.exists(tree_file):
                        tree_files.append(tree_file)

                        if i % 10 == 0:
                            logger.info(f"Processed {i}/{len(window_files)} windows...")
                    else:
                        logger.warning(f"Tree file not found for {window_file}")
                        failed_windows.append(window_file)
                else:
                    logger.warning(f"IQ-TREE failed for {window_file}: {result.stderr}")
                    failed_windows.append(window_file)

            except subprocess.TimeoutExpired:
                logger.warning(f"Timeout for window {window_file}")
                failed_windows.append(window_file)

        logger.info(
            f"Successfully inferred {len(tree_files)} trees ({len(failed_windows)} failed)"
        )

        if failed_windows:
            logger.info(f"Failed windows: {len(failed_windows)}")

        return tree_files

    def validate_rooting_consistency(self, tree_files: List[str]) -> Dict:
        """
        Validate rooting consistency using the external validator script.

        Args:
            tree_files: List of tree files to validate

        Returns:
            Dictionary with validation results
        """
        logger.info("Validating rooting consistency...")

        # Use the rooting consistency validator
        validator_script = Path(__file__).parent / "validate_rooting_consistency.py"

        if not validator_script.exists():
            logger.warning(
                "Rooting consistency validator not found. Performing basic validation."
            )
            return self._basic_rooting_validation(tree_files)

        # Run the validator
        validation_dir = self.results_dir / "rooting_validation"
        validation_dir.mkdir(exist_ok=True)

        try:
            cmd = [
                "python3",
                str(validator_script),
                "--output-dir",
                str(validation_dir),
            ] + tree_files

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            if result.returncode == 0:
                # Load results
                metrics_file = validation_dir / "consistency_metrics.json"
                if metrics_file.exists():
                    with open(metrics_file, "r") as f:
                        return json.load(f)
            else:
                logger.warning(f"Validation script failed: {result.stderr}")

        except subprocess.TimeoutExpired:
            logger.warning("Validation script timed out")
        except Exception as e:
            logger.warning(f"Validation script error: {e}")

        # Fallback to basic validation
        return self._basic_rooting_validation(tree_files)

    def _basic_rooting_validation(self, tree_files: List[str]) -> Dict:
        """Basic rooting validation without external dependencies."""
        from Bio import Phylo

        rooted_count = 0
        total_count = 0

        for tree_file in tree_files:
            try:
                tree = Phylo.read(tree_file, "newick")
                total_count += 1
                if tree.rooted:
                    rooted_count += 1
            except Exception:
                pass

        if total_count > 0:
            consistency_score = rooted_count / total_count
        else:
            consistency_score = 0

        return {
            "overall_consistency_score": consistency_score,
            "num_trees": total_count,
            "num_rooted_trees": rooted_count,
            "rooted_percentage": (rooted_count / total_count * 100)
            if total_count > 0
            else 0,
        }

    def run_complete_analysis(
        self, alignment_file: str, window_size: int = 400, step_size: int = 5
    ) -> Dict:
        """
        Run the complete constraint-based rooting analysis.

        Args:
            alignment_file: Path to the input alignment
            window_size: Window size in base pairs
            step_size: Step size in base pairs

        Returns:
            Dictionary with complete analysis results
        """
        logger.info("Starting complete constraint-based rooting analysis...")

        try:
            # Step 1: Create backbone constraint tree
            constraint_tree = self.create_backbone_constraint_tree(alignment_file)

            # Step 2: Generate sliding windows
            window_files = self.generate_sliding_windows(
                alignment_file, window_size, step_size
            )

            if not window_files:
                raise RuntimeError("No valid sliding windows generated")

            # Step 3: Infer constraint-guided trees
            tree_files = self.infer_constraint_guided_trees(
                window_files, constraint_tree
            )

            if not tree_files:
                raise RuntimeError("No trees successfully inferred")

            # Step 4: Validate rooting consistency
            validation_results = self.validate_rooting_consistency(tree_files)

            # Step 5: Combine all trees into a single file for easy access
            combined_tree_file = self.results_dir / "all_constraint_trees.newick"
            with open(combined_tree_file, "w") as out_f:
                for tree_file in tree_files:
                    with open(tree_file, "r") as in_f:
                        out_f.write(in_f.read().strip() + "\n")

            # Compile results
            results = {
                "analysis_type": "constraint_based_rooting",
                "input_alignment": alignment_file,
                "window_size": window_size,
                "step_size": step_size,
                "backbone_constraint_tree": constraint_tree,
                "num_windows_generated": len(window_files),
                "num_trees_inferred": len(tree_files),
                "combined_tree_file": str(combined_tree_file),
                "validation_results": validation_results,
                "success_rate": len(tree_files) / len(window_files) * 100
                if window_files
                else 0,
            }

            # Save results
            results_file = self.results_dir / "constraint_rooting_results.json"
            with open(results_file, "w") as f:
                json.dump(results, f, indent=2, default=str)

            # Generate summary report
            self._generate_summary_report(results)

            logger.info("Complete analysis finished successfully!")
            return results

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise

    def _generate_summary_report(self, results: Dict):
        """Generate a summary report of the analysis."""
        report_file = self.results_dir / "CONSTRAINT_ROOTING_SUMMARY.md"

        with open(report_file, "w") as f:
            f.write("# Constraint-Based Rooting Analysis Summary\n\n")
            f.write(f"**Analysis Date**: {pd.Timestamp.now()}\n")
            f.write(f"**Method**: IQ-TREE2 with backbone constraint trees\n")
            f.write(f"**Input Alignment**: {results['input_alignment']}\n\n")

            f.write("## Parameters\n\n")
            f.write(f"- **Window Size**: {results['window_size']} bp\n")
            f.write(f"- **Step Size**: {results['step_size']} bp\n")
            f.write(
                f"- **Overlap**: {((results['window_size'] - results['step_size']) / results['window_size'] * 100):.1f}%\n\n"
            )

            f.write("## Results\n\n")
            f.write(f"- **Windows Generated**: {results['num_windows_generated']}\n")
            f.write(
                f"- **Trees Successfully Inferred**: {results['num_trees_inferred']}\n"
            )
            f.write(f"- **Success Rate**: {results['success_rate']:.1f}%\n\n")

            # Validation results
            validation = results.get("validation_results", {})
            if validation:
                consistency_score = validation.get("overall_consistency_score", 0)
                f.write(f"- **Rooting Consistency Score**: {consistency_score:.3f}\n")
                f.write(
                    f"- **Rooted Trees**: {validation.get('num_rooted_trees', 0)}/{validation.get('num_trees', 0)}\n\n"
                )

                # Interpretation
                f.write("## Interpretation\n\n")
                if consistency_score >= 0.95:
                    f.write(
                        "🎉 **Excellent Results!** Constraint-based rooting achieved >95% consistency.\n"
                    )
                    f.write("✅ This approach is working very well for your dataset.\n")
                    f.write("✅ Results are suitable for publication.\n\n")
                elif consistency_score >= 0.85:
                    f.write(
                        "✅ **Good Results!** Constraint-based rooting achieved good consistency (85-95%).\n"
                    )
                    f.write(
                        "🔧 Minor optimizations could improve consistency further.\n\n"
                    )
                elif consistency_score >= 0.70:
                    f.write(
                        "⚠️ **Moderate Results.** Some improvement in consistency achieved.\n"
                    )
                    f.write(
                        "🔧 Consider adjusting window parameters or trying non-reversible models.\n\n"
                    )
                else:
                    f.write(
                        "❌ **Limited Improvement.** Constraint approach didn't significantly improve consistency.\n"
                    )
                    f.write(
                        "🔧 **Recommendation**: Try direct rooted inference with non-reversible models.\n\n"
                    )

            f.write("## Files Generated\n\n")
            f.write(
                f"- **Backbone Constraint Tree**: `{Path(results['backbone_constraint_tree']).name}`\n"
            )
            f.write(
                f"- **Combined Trees**: `{Path(results['combined_tree_file']).name}`\n"
            )
            f.write(
                "- **Individual Tree Files**: `sliding_windows/window_*_constrained.treefile`\n"
            )
            f.write(
                "- **Validation Report**: `rooting_validation/rooting_consistency_report.md`\n\n"
            )

            f.write("## Next Steps\n\n")
            if validation.get("overall_consistency_score", 0) >= 0.85:
                f.write("1. ✅ Use these trees for your sliding window analysis\n")
                f.write("2. 📊 Compare topological congruence with previous results\n")
                f.write("3. 📝 Document this approach for your publication\n")
            else:
                f.write(
                    "1. 🔧 Try the non-reversible model approach (direct rooted inference)\n"
                )
                f.write("2. 🔧 Experiment with different window sizes (300-500 bp)\n")
                f.write(
                    "3. 🔧 Consider professional tools like RDP4/RDP5 for publication quality\n"
                )

            f.write("\n---\n")
            f.write("*Generated by Constraint Rooting Implementation v1.0*\n")

        logger.info(f"Summary report generated: {report_file}")


def main():
    """Main execution function."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Implement constraint-based rooting for sliding window analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with default parameters
    python implement_constraint_rooting.py --alignment ingroup_only_alignment.fasta
    
    # Custom window parameters
    python implement_constraint_rooting.py --alignment ingroup_only_alignment.fasta --window-size 500 --step-size 10
        """,
    )

    parser.add_argument(
        "--alignment", required=True, help="Input alignment file (FASTA format)"
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=400,
        help="Window size in base pairs (default: 400)",
    )
    parser.add_argument(
        "--step-size", type=int, default=5, help="Step size in base pairs (default: 5)"
    )
    parser.add_argument(
        "--output-dir",
        default="results_constraint_rooting",
        help="Output directory (default: results_constraint_rooting)",
    )

    args = parser.parse_args()

    # Validate input file
    if not os.path.exists(args.alignment):
        print(f"Error: Alignment file not found: {args.alignment}")
        sys.exit(1)

    print("🧬 Constraint-Based Rooting Implementation for Sliding Window Analysis")
    print(f"📁 Input: {args.alignment}")
    print(f"📊 Window: {args.window_size} bp, Step: {args.step_size} bp")
    print(f"📂 Output: {args.output_dir}")
    print()

    try:
        # Initialize and run analysis
        implementation = ConstraintRootingImplementation(args.output_dir)
        results = implementation.run_complete_analysis(
            args.alignment, args.window_size, args.step_size
        )

        # Display results
        print("🎉 Analysis completed successfully!")
        print()
        print("📊 Summary:")
        print(f"   Windows Generated: {results['num_windows_generated']}")
        print(f"   Trees Inferred: {results['num_trees_inferred']}")
        print(f"   Success Rate: {results['success_rate']:.1f}%")

        validation = results.get("validation_results", {})
        if validation:
            consistency_score = validation.get("overall_consistency_score", 0)
            print(f"   Rooting Consistency: {consistency_score:.3f}")

            if consistency_score >= 0.95:
                print("   🎉 Excellent consistency achieved!")
            elif consistency_score >= 0.85:
                print("   ✅ Good consistency achieved!")
            elif consistency_score >= 0.70:
                print("   ⚠️ Moderate consistency achieved")
            else:
                print("   ❌ Limited consistency improvement")

        print()
        print(f"📁 Results saved to: {args.output_dir}/")
        print(f"📋 Summary report: {args.output_dir}/CONSTRAINT_ROOTING_SUMMARY.md")
        print(f"🌳 Combined trees: {args.output_dir}/all_constraint_trees.newick")

    except KeyboardInterrupt:
        print("\n❌ Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
