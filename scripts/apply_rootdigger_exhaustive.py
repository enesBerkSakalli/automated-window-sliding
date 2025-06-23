#!/usr/bin/env python3
"""
Apply RootDigger Exhaustive Search to Sliding Window Trees

This script applies RootDigger's exhaustive search to find optimal root positions
for all sliding window trees, ensuring consistent and statistically optimal rooting.

Usage:
    python apply_rootdigger_exhaustive.py --trees results_constraint_test/sliding_windows/window_*_constrained.treefile
    python apply_rootdigger_exhaustive.py --alignment ingroup_only_alignment.fasta --auto-find-trees

Author: Enhanced Norovirus Pipeline
Date: January 2025
"""

import os
import sys
import subprocess
import argparse
import logging
from pathlib import Path
from typing import List, Dict
import glob
import json
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("rootdigger_exhaustive.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


class RootDiggerExhaustive:
    """
    Apply RootDigger exhaustive search to sliding window trees.
    """

    def __init__(self, output_dir: str = "results_rootdigger_exhaustive"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Check RootDigger availability
        if not self._check_rootdigger():
            raise RuntimeError(
                "RootDigger is required but not found. Please install RootDigger."
            )

        logger.info(
            f"RootDigger exhaustive rooting initialized. Results: {self.output_dir}"
        )

    def _check_rootdigger(self) -> bool:
        """Check if RootDigger is available."""
        try:
            result = subprocess.run(
                ["rootdigger", "--help"], capture_output=True, text=True
            )
            return result.returncode == 0
        except FileNotFoundError:
            return False

    def find_sliding_window_files(self, pattern: str = None) -> tuple:
        """
        Find sliding window alignment and tree files.

        Args:
            pattern: Custom pattern for finding files

        Returns:
            Tuple of (alignment_files, tree_files)
        """
        if pattern:
            tree_files = glob.glob(pattern)
        else:
            # Look for common patterns
            patterns = [
                "results_constraint_test/sliding_windows/window_*_constrained.treefile",
                "results_*/sliding_windows/window_*.treefile",
                "sliding_windows/window_*.treefile",
                "window_*.treefile",
            ]

            tree_files = []
            for pat in patterns:
                found = glob.glob(pat)
                if found:
                    tree_files.extend(found)
                    logger.info(f"Found {len(found)} trees with pattern: {pat}")
                    break

        # Find corresponding alignment files
        alignment_files = []
        for tree_file in tree_files:
            # Try to find corresponding alignment file
            tree_path = Path(tree_file)
            base_name = tree_path.stem.replace("_constrained", "").replace(
                "_rooted", ""
            )

            # Look for alignment with same base name
            alignment_patterns = [
                tree_path.parent / f"{base_name}.fasta",
                tree_path.parent / f"{base_name}.fa",
                tree_path.parent / f"{base_name}.phy",
            ]

            alignment_file = None
            for align_pat in alignment_patterns:
                if align_pat.exists():
                    alignment_file = str(align_pat)
                    break

            if alignment_file:
                alignment_files.append(alignment_file)
            else:
                logger.warning(f"No alignment found for tree: {tree_file}")

        logger.info(
            f"Found {len(tree_files)} trees and {len(alignment_files)} alignments"
        )
        return alignment_files, tree_files

    def apply_rootdigger_to_tree(
        self, alignment_file: str, tree_file: str, threads: int = 1
    ) -> str:
        """
        Apply RootDigger exhaustive search to a single tree.

        Args:
            alignment_file: Path to alignment file
            tree_file: Path to tree file
            threads: Number of threads to use

        Returns:
            Path to rooted tree file
        """
        tree_name = Path(tree_file).stem
        output_prefix = self.output_dir / f"{tree_name}_rootdigger"

        cmd = [
            "rootdigger",
            "--msa",
            alignment_file,
            "--tree",
            tree_file,
            "--exhaustive",  # Use exhaustive search
            "--threads",
            str(threads),
            "--prefix",
            str(output_prefix),
            "--silent",  # Reduce output verbosity
        ]

        logger.debug(f"Running RootDigger: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=1800
            )  # 30-minute timeout

            if result.returncode == 0:
                # RootDigger outputs the rooted tree with .rooted.tree extension
                rooted_tree = f"{output_prefix}.rooted.tree"
                if os.path.exists(rooted_tree):
                    return rooted_tree
                else:
                    logger.warning(f"Expected rooted tree not found: {rooted_tree}")
                    return None
            else:
                logger.error(f"RootDigger failed for {tree_file}: {result.stderr}")
                return None

        except subprocess.TimeoutExpired:
            logger.warning(f"RootDigger timed out for {tree_file}")
            return None
        except Exception as e:
            logger.error(f"RootDigger error for {tree_file}: {e}")
            return None

    def process_all_trees(
        self,
        alignment_files: List[str],
        tree_files: List[str],
        threads: int = 1,
        max_parallel: int = 1,
    ) -> List[str]:
        """
        Process all sliding window trees with RootDigger.

        Args:
            alignment_files: List of alignment files
            tree_files: List of tree files
            threads: Threads per RootDigger process
            max_parallel: Maximum parallel RootDigger processes

        Returns:
            List of successfully rooted tree files
        """
        logger.info(
            f"Processing {len(tree_files)} trees with RootDigger exhaustive search..."
        )

        if len(alignment_files) != len(tree_files):
            logger.warning("Mismatch between number of alignments and trees")
            # Pad with None for missing alignments
            while len(alignment_files) < len(tree_files):
                alignment_files.append(None)

        rooted_trees = []
        failed_trees = []

        for i, (alignment_file, tree_file) in enumerate(
            zip(alignment_files, tree_files), 1
        ):
            if alignment_file is None:
                logger.warning(f"No alignment for tree {tree_file}, skipping")
                failed_trees.append(tree_file)
                continue

            logger.info(
                f"Processing tree {i}/{len(tree_files)}: {Path(tree_file).name}"
            )

            rooted_tree = self.apply_rootdigger_to_tree(
                alignment_file, tree_file, threads
            )

            if rooted_tree:
                rooted_trees.append(rooted_tree)
            else:
                failed_trees.append(tree_file)

            # Progress update
            if i % 5 == 0:
                success_rate = len(rooted_trees) / i * 100
                logger.info(
                    f"Progress: {i}/{len(tree_files)} ({success_rate:.1f}% success rate)"
                )

        logger.info(
            f"RootDigger processing completed: {len(rooted_trees)} successful, {len(failed_trees)} failed"
        )

        if failed_trees:
            logger.info(f"Failed trees: {[Path(f).name for f in failed_trees]}")

        return rooted_trees

    def create_combined_tree_file(self, rooted_trees: List[str]) -> str:
        """Create a combined file with all rooted trees."""
        combined_file = self.output_dir / "all_rootdigger_trees.newick"

        with open(combined_file, "w") as out_f:
            for tree_file in rooted_trees:
                with open(tree_file, "r") as in_f:
                    tree_content = in_f.read().strip()
                    out_f.write(tree_content + "\n")

        logger.info(f"Combined rooted trees saved to: {combined_file}")
        return str(combined_file)

    def validate_rooting_consistency(self, rooted_trees: List[str]) -> Dict:
        """Validate rooting consistency using the external validator."""
        logger.info("Validating RootDigger rooting consistency...")

        validator_script = Path(__file__).parent / "validate_rooting_consistency.py"

        if not validator_script.exists():
            logger.warning("Rooting consistency validator not found")
            return {}

        validation_dir = self.output_dir / "rooting_validation"
        validation_dir.mkdir(exist_ok=True)

        try:
            cmd = [
                "python3",
                str(validator_script),
                "--output-dir",
                str(validation_dir),
            ] + rooted_trees

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            if result.returncode == 0:
                # Load results
                metrics_file = validation_dir / "consistency_metrics.json"
                if metrics_file.exists():
                    with open(metrics_file, "r") as f:
                        return json.load(f)
            else:
                logger.warning(f"Validation failed: {result.stderr}")

        except Exception as e:
            logger.warning(f"Validation error: {e}")

        return {}

    def generate_summary_report(self, results: Dict):
        """Generate summary report of RootDigger analysis."""
        report_file = self.output_dir / "ROOTDIGGER_EXHAUSTIVE_SUMMARY.md"

        with open(report_file, "w") as f:
            f.write("# RootDigger Exhaustive Search Summary\n\n")
            f.write(f"**Analysis Date**: {pd.Timestamp.now()}\n")
            f.write(
                "**Method**: RootDigger exhaustive search for optimal root placement\n"
            )
            f.write(f"**RootDigger Version**: v1.8.0+\n\n")

            f.write("## Results\n\n")
            f.write(f"- **Input Trees**: {results.get('num_input_trees', 'N/A')}\n")
            f.write(
                f"- **Successfully Rooted**: {results.get('num_rooted_trees', 'N/A')}\n"
            )
            f.write(f"- **Success Rate**: {results.get('success_rate', 0):.1f}%\n\n")

            # Validation results
            validation = results.get("validation_results", {})
            if validation:
                consistency_score = validation.get("overall_consistency_score", 0)
                f.write(f"- **Rooting Consistency Score**: {consistency_score:.3f}\n")
                f.write(
                    f"- **Root Depth Consistency**: {validation.get('root_depth_consistency', 'N/A'):.3f}\n"
                )
                f.write(
                    f"- **Split Pattern Consistency**: {validation.get('split_pattern_consistency', 'N/A'):.3f}\n\n"
                )

                # Interpretation
                f.write("## Interpretation\n\n")
                if consistency_score >= 0.95:
                    f.write(
                        "🎉 **Outstanding Results!** RootDigger achieved >95% rooting consistency.\n"
                    )
                    f.write(
                        "✅ Exhaustive search found highly consistent optimal root positions.\n"
                    )
                    f.write("✅ Results are publication-quality.\n\n")
                elif consistency_score >= 0.85:
                    f.write(
                        "✅ **Excellent Results!** RootDigger achieved high consistency (85-95%).\n"
                    )
                    f.write(
                        "🔬 Exhaustive search significantly improved root placement.\n\n"
                    )
                elif consistency_score >= 0.70:
                    f.write(
                        "⚠️ **Good Results.** Moderate improvement in rooting consistency.\n"
                    )
                    f.write(
                        "🔧 Consider combining with constraint trees for better results.\n\n"
                    )
                else:
                    f.write(
                        "❌ **Limited Improvement.** Root positions still highly variable.\n"
                    )
                    f.write(
                        "🔧 **Recommendation**: Check alignment quality and window parameters.\n\n"
                    )

            f.write("## RootDigger Method\n\n")
            f.write("RootDigger uses exhaustive search to:\n")
            f.write("1. Test all possible root positions on each tree\n")
            f.write("2. Optimize likelihood for each root placement\n")
            f.write("3. Select the statistically optimal root position\n")
            f.write("4. Provide rooted trees with maximum likelihood support\n\n")

            f.write("## Files Generated\n\n")
            f.write(f"- **Combined Rooted Trees**: `all_rootdigger_trees.newick`\n")
            f.write("- **Individual Rooted Trees**: `window_*_rootdigger.rooted`\n")
            f.write(
                "- **Validation Report**: `rooting_validation/rooting_consistency_report.md`\n"
            )
            f.write(
                "- **Detailed Metrics**: `rooting_validation/consistency_metrics.json`\n\n"
            )

            f.write("## Next Steps\n\n")
            if validation.get("overall_consistency_score", 0) >= 0.85:
                f.write(
                    "1. ✅ Use these optimally rooted trees for downstream analysis\n"
                )
                f.write("2. 📊 Compare topological congruence with previous methods\n")
                f.write("3. 📝 Document RootDigger methodology for publication\n")
                f.write("4. 🔬 Consider these as your final rooted trees\n")
            else:
                f.write(
                    "1. 🔧 Combine RootDigger with constraint trees for hybrid approach\n"
                )
                f.write(
                    "2. 🔧 Evaluate alignment quality and remove problematic sequences\n"
                )
                f.write("3. 🔧 Consider larger window sizes or different step sizes\n")
                f.write(
                    "4. 🔬 Professional tools (RDP4/RDP5) may still be beneficial\n"
                )

            f.write("\n---\n")
            f.write("*Generated by RootDigger Exhaustive Analysis v1.0*\n")

        logger.info(f"Summary report generated: {report_file}")

    def run_complete_analysis(
        self,
        alignment_files: List[str] = None,
        tree_files: List[str] = None,
        pattern: str = None,
        threads: int = 1,
    ) -> Dict:
        """
        Run complete RootDigger exhaustive analysis.

        Args:
            alignment_files: List of alignment files (optional, will auto-find)
            tree_files: List of tree files (optional, will auto-find)
            pattern: Pattern for finding tree files
            threads: Number of threads per RootDigger process

        Returns:
            Dictionary with analysis results
        """
        logger.info("Starting RootDigger exhaustive rooting analysis...")

        try:
            # Find files if not provided
            if not tree_files:
                alignment_files, tree_files = self.find_sliding_window_files(pattern)

            if not tree_files:
                raise RuntimeError("No tree files found for processing")

            # Process all trees with RootDigger
            rooted_trees = self.process_all_trees(alignment_files, tree_files, threads)

            if not rooted_trees:
                raise RuntimeError("No trees successfully rooted by RootDigger")

            # Create combined tree file
            combined_file = self.create_combined_tree_file(rooted_trees)

            # Validate rooting consistency
            validation_results = self.validate_rooting_consistency(rooted_trees)

            # Compile results
            results = {
                "method": "rootdigger_exhaustive",
                "num_input_trees": len(tree_files),
                "num_rooted_trees": len(rooted_trees),
                "success_rate": len(rooted_trees) / len(tree_files) * 100
                if tree_files
                else 0,
                "combined_tree_file": combined_file,
                "rooted_tree_files": rooted_trees,
                "validation_results": validation_results,
            }

            # Save results
            results_file = self.output_dir / "rootdigger_results.json"
            with open(results_file, "w") as f:
                json.dump(results, f, indent=2, default=str)

            # Generate summary report
            self.generate_summary_report(results)

            logger.info("RootDigger exhaustive analysis completed successfully!")
            return results

        except Exception as e:
            logger.error(f"RootDigger analysis failed: {e}")
            raise


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Apply RootDigger exhaustive search to sliding window trees",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Auto-find trees from constraint rooting results
    python apply_rootdigger_exhaustive.py --auto-find
    
    # Specify tree pattern
    python apply_rootdigger_exhaustive.py --pattern "results_*/sliding_windows/window_*_constrained.treefile"
    
    # Specify individual files
    python apply_rootdigger_exhaustive.py --trees tree1.treefile tree2.treefile --alignments align1.fasta align2.fasta
        """,
    )

    parser.add_argument("--trees", nargs="+", help="Tree files to process")
    parser.add_argument("--alignments", nargs="+", help="Corresponding alignment files")
    parser.add_argument(
        "--pattern",
        help='Pattern for finding tree files (e.g., "results_*/window_*.treefile")',
    )
    parser.add_argument(
        "--auto-find", action="store_true", help="Auto-find sliding window files"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads per RootDigger process (default: 1)",
    )
    parser.add_argument(
        "--output-dir", default="results_rootdigger_exhaustive", help="Output directory"
    )

    args = parser.parse_args()

    # Validate arguments
    if not (args.trees or args.pattern or args.auto_find):
        print("Error: Must specify --trees, --pattern, or --auto-find")
        parser.print_help()
        sys.exit(1)

    print("🔍 RootDigger Exhaustive Search for Sliding Window Trees")
    print("=" * 60)
    print("🚀 Finding optimal root positions using exhaustive likelihood search")
    print()

    try:
        # Initialize RootDigger processor
        rootdigger = RootDiggerExhaustive(args.output_dir)

        # Run analysis
        if args.auto_find:
            results = rootdigger.run_complete_analysis(threads=args.threads)
        elif args.pattern:
            results = rootdigger.run_complete_analysis(
                pattern=args.pattern, threads=args.threads
            )
        else:
            results = rootdigger.run_complete_analysis(
                alignment_files=args.alignments,
                tree_files=args.trees,
                threads=args.threads,
            )

        # Display results
        print("🎉 RootDigger analysis completed!")
        print()
        print("📊 Summary:")
        print(f"   Input Trees: {results['num_input_trees']}")
        print(f"   Successfully Rooted: {results['num_rooted_trees']}")
        print(f"   Success Rate: {results['success_rate']:.1f}%")

        validation = results.get("validation_results", {})
        if validation:
            consistency_score = validation.get("overall_consistency_score", 0)
            print(f"   Rooting Consistency: {consistency_score:.3f}")

            if consistency_score >= 0.95:
                print("   🎉 Outstanding consistency achieved!")
            elif consistency_score >= 0.85:
                print("   ✅ Excellent consistency achieved!")
            elif consistency_score >= 0.70:
                print("   ⚠️ Good consistency achieved")
            else:
                print("   ❌ Limited consistency improvement")

        print()
        print(f"📁 Results: {args.output_dir}/")
        print(f"📋 Report: {args.output_dir}/ROOTDIGGER_EXHAUSTIVE_SUMMARY.md")
        print(f"🌳 Trees: {args.output_dir}/all_rootdigger_trees.newick")

    except KeyboardInterrupt:
        print("\n❌ Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
