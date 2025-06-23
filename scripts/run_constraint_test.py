#!/usr/bin/env python3
"""
Quick Run: Test Constraint-Based Rooting on Norovirus Dataset

This script runs the constraint-based rooting implementation on the current
norovirus dataset and compares results with previous approaches.

Usage:
    python run_constraint_test.py
"""

import subprocess
import sys
from pathlib import Path
import json


def main():
    """Run constraint-based rooting test."""

    # Check if alignment file exists
    alignment_file = "ingroup_only_alignment.fasta"
    if not Path(alignment_file).exists():
        print(f"❌ Alignment file not found: {alignment_file}")
        print(
            "Please ensure you're in the correct directory with the norovirus alignment."
        )
        return False

    print("🧬 Testing Constraint-Based Rooting on Norovirus Dataset")
    print("=" * 60)

    # Run the constraint rooting implementation
    print("🚀 Starting constraint-based analysis...")
    print("   This may take 10-30 minutes depending on your system.")
    print()

    try:
        cmd = [
            "python3",
            "implement_constraint_rooting.py",
            "--alignment",
            alignment_file,
            "--window-size",
            "400",  # Optimized from previous analysis
            "--step-size",
            "5",  # High overlap for consistency
            "--output-dir",
            "results_constraint_test",
        ]

        print(f"Running: {' '.join(cmd)}")
        print()

        result = subprocess.run(cmd, check=True)

        if result.returncode == 0:
            print()
            print("✅ Constraint-based analysis completed!")

            # Load and display results
            results_file = Path(
                "results_constraint_test/constraint_rooting_results.json"
            )
            if results_file.exists():
                with open(results_file, "r") as f:
                    results = json.load(f)

                print()
                print("📊 Quick Results Summary:")
                print(f"   Windows: {results.get('num_windows_generated', 'N/A')}")
                print(f"   Trees: {results.get('num_trees_inferred', 'N/A')}")
                print(f"   Success Rate: {results.get('success_rate', 0):.1f}%")

                validation = results.get("validation_results", {})
                if validation:
                    consistency = validation.get("overall_consistency_score", 0)
                    print(f"   Consistency: {consistency:.3f}")

                    if consistency >= 0.95:
                        print("   🎉 EXCELLENT - Major improvement achieved!")
                    elif consistency >= 0.85:
                        print("   ✅ GOOD - Significant improvement!")
                    elif consistency >= 0.70:
                        print("   ⚠️ MODERATE - Some improvement")
                    else:
                        print("   ❌ LIMITED - Try alternative approach")

            print()
            print("📁 Results saved to: results_constraint_test/")
            print(
                "📋 Full report: results_constraint_test/CONSTRAINT_ROOTING_SUMMARY.md"
            )
            print("🌳 Trees: results_constraint_test/all_constraint_trees.newick")

            # Compare with previous results if available
            print()
            print("🔍 Comparison with Previous Results:")

            prev_rf = 0.87  # From previous improved strategy
            if validation and "overall_consistency_score" in validation:
                new_consistency = validation["overall_consistency_score"]

                if new_consistency > 0.95:
                    print(
                        f"   📈 Rooting consistency: {new_consistency:.3f} (target: >0.95) ✅"
                    )
                elif new_consistency > prev_rf:
                    improvement = ((new_consistency - prev_rf) / prev_rf) * 100
                    print(f"   📈 Improvement over previous: +{improvement:.1f}%")
                else:
                    print(
                        f"   📊 Consistency: {new_consistency:.3f} (previous RF: {prev_rf})"
                    )

            return True

    except subprocess.CalledProcessError as e:
        print(f"❌ Analysis failed with return code {e.returncode}")
        return False
    except FileNotFoundError:
        print("❌ implement_constraint_rooting.py not found")
        print("Please ensure the script is in the current directory.")
        return False
    except KeyboardInterrupt:
        print("\n❌ Analysis interrupted by user")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False


if __name__ == "__main__":
    success = main()

    if success:
        print()
        print("🎯 Next Steps:")
        print("1. Review the summary report for detailed analysis")
        print("2. If consistency >0.95, use these trees for your research")
        print("3. If consistency <0.95, try the non-reversible model approach")
        print("4. Compare topological congruence with previous results")
        print()
        print("📚 For publication-quality analysis, consider:")
        print("   - RDP4/RDP5 for professional recombination analysis")
        print("   - SlidingBayes for Bayesian sliding window analysis")
        sys.exit(0)
    else:
        print()
        print("🔧 Troubleshooting:")
        print("1. Ensure IQ-TREE2 is installed and in PATH")
        print("2. Check that input alignment file exists")
        print("3. Verify sufficient disk space and memory")
        sys.exit(1)
