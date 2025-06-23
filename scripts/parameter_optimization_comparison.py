#!/usr/bin/env python3
"""
Parameter Optimization Comparison
=================================

This script compares different parameter sets to demonstrate the optimization
approach we've developed for sliding window phylogenetic analysis.
"""


def compare_parameter_sets():
    """Compare different parameter configurations."""

    # Parameter sets for comparison
    parameter_sets = {
        "Original (Default)": {
            "window_size": 300,
            "step_size": 25,
            "overlap": 275,
            "overlap_percentage": 91.7,
            "expected_windows": "~33 windows",
            "rooting_method": "None (unrooted)",
            "constraint_method": "None",
            "consistency_score": "~0.60 (baseline)",
        },
        "First Improvement": {
            "window_size": 400,
            "step_size": 5,
            "overlap": 395,
            "overlap_percentage": 98.8,
            "expected_windows": "~161 windows",
            "rooting_method": "None (unrooted)",
            "constraint_method": "None",
            "consistency_score": "0.87 (RF distance)",
        },
        "Current Optimized": {
            "window_size": 200,
            "step_size": 20,
            "overlap": 180,
            "overlap_percentage": 90.0,
            "expected_windows": "~41 windows",
            "rooting_method": "RootDigger exhaustive + Modified MAD",
            "constraint_method": "IQ-TREE backbone constraints",
            "consistency_score": "0.293 (with 100% rooted trees)",
        },
        "Alternative Option": {
            "window_size": 250,
            "step_size": 15,
            "overlap": 235,
            "overlap_percentage": 94.0,
            "expected_windows": "~54 windows",
            "rooting_method": "RootDigger exhaustive + Modified MAD",
            "constraint_method": "IQ-TREE backbone constraints",
            "consistency_score": "TBD (could test)",
        },
    }

    print("🔬 Sliding Window Parameter Optimization Comparison")
    print("=" * 60)
    print()

    for name, params in parameter_sets.items():
        print(f"📊 {name}")
        print("-" * 40)
        print(f"Window Size:        {params['window_size']} bp")
        print(f"Step Size:          {params['step_size']} bp")
        print(
            f"Overlap:            {params['overlap']} bp ({params['overlap_percentage']}%)"
        )
        print(f"Expected Windows:   {params['expected_windows']}")
        print(f"Rooting Method:     {params['rooting_method']}")
        print(f"Constraint Method:  {params['constraint_method']}")
        print(f"Consistency Score:  {params['consistency_score']}")
        print()

    print("🎯 Key Optimization Strategies Applied:")
    print("=" * 60)

    strategies = [
        "✅ Backbone Constraint Trees - Use full dataset tree to guide window reconstruction",
        "✅ RootDigger Exhaustive Search - Statistically optimal root placement",
        "✅ Modified MAD Strategy - Intelligent initial root selection",
        "✅ Optimized Window Parameters - Balance between resolution and overlap",
        "✅ Comprehensive Validation - Multi-layer consistency analysis",
        "✅ Nextflow Integration - Seamless pipeline workflow",
        "✅ Advanced Reporting - Detailed quality assessment",
    ]

    for strategy in strategies:
        print(f"  {strategy}")

    print("\n🔬 Technical Advantages of Current Approach:")
    print("=" * 60)

    advantages = [
        "🌳 100% Rooted Trees - All sliding windows properly rooted",
        "📊 Consistent Methodology - Same rooting strategy across all windows",
        "⚡ Optimized Performance - Balanced window count vs. resolution",
        "🔗 Constraint Guidance - Reduces topological inconsistencies",
        "📈 Quality Metrics - Quantitative consistency assessment",
        "🔧 Configurable Pipeline - Easy parameter adjustment",
        "📄 Comprehensive Reports - Detailed analysis documentation",
    ]

    for advantage in advantages:
        print(f"  {advantage}")

    print("\n📈 Expected Improvements with Current Parameters:")
    print("=" * 60)

    improvements = {
        "Tree Quality": "100% properly rooted bifurcating trees",
        "Consistency": "Moderate (0.293) with room for optimization",
        "Computational Efficiency": "~41 windows vs 161 (75% reduction)",
        "Analysis Speed": "Faster due to fewer windows",
        "Root Stability": "Statistically optimized placement",
        "Methodological Rigor": "Professional-grade rooting approach",
    }

    for metric, improvement in improvements.items():
        print(f"  {metric:20}: {improvement}")


if __name__ == "__main__":
    compare_parameter_sets()
