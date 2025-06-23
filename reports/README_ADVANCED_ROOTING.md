# Advanced Rooting Tools for Sliding Window Phylogenetic Analysis

## 🎯 Executive Summary

This document provides a comprehensive guide to state-of-the-art tools and strategies for achieving consistent rooting in sliding window phylogenetic analyses, based on extensive research of the 2024-2025 software landscape.

## 🔬 Research Summary

Our investigation identified several advanced approaches for maintaining consistent rooting across sliding window trees while preserving bifurcating topology and maximizing topological congruence:

### ⭐ Top-Tier Solutions

1. **IQ-TREE 2.3+ with Constraint Trees** - Best current option
2. **RAxML-NG with Topology Constraints** - Excellent alternative
3. **Non-reversible Models (UNREST/NONREV)** - Direct rooted inference
4. **Statistical Root Testing** - Confidence assessment
5. **Professional Tools** - RDP4/RDP5, SlidingBayes

## 🚀 Quick Start - Immediate Implementation

### Option 1: Test Constraint-Based Rooting (Recommended)

```bash
# Quick test on your norovirus dataset
python run_constraint_test.py
```

This will:
- ✅ Create a backbone constraint tree from your full alignment
- ✅ Apply constraints to all sliding window trees
- ✅ Validate rooting consistency
- ✅ Generate a comprehensive report

### Option 2: Full Advanced Pipeline

```bash
# Complete analysis with multiple strategies
python advanced_rooting_pipeline.py --alignment ingroup_only_alignment.fasta --strategy constraint_guided
```

### Option 3: Manual IQ-TREE Implementation

```bash
# Step 1: Create backbone constraint tree
iqtree2 -s ingroup_only_alignment.fasta -m MFP -B 1000 --prefix backbone_constraint

# Step 2: Apply to sliding windows (example for one window)
iqtree2 -s window_001.fasta -g backbone_constraint.treefile -B 1000 --prefix window_001_constrained
```

## 🛠️ Advanced Tools and Methods

### IQ-TREE 2.3+ Features

**Constraint Trees** (`-g constraint.treefile`)
- Force consistent topology across windows
- Maintain root position relationships
- High statistical support with bootstrap

**Non-reversible Models** (`--model-joint UNREST`)
- Direct inference of rooted trees
- No post-hoc rooting required
- Statistical root confidence (rootstrap)

**Root Testing** (`--root-test`)
- Statistical validation of root positions
- AU, KH, SH tests for significance
- Quantitative root confidence scores

### RAxML-NG Integration

```bash
# Constraint-based approach
raxml-ng --msa window.fasta --tree-constraint backbone.tree --model GTR+G --prefix constrained_window
```

### Professional-Grade Solutions

**RDP4/RDP5** (Recombination Detection Program)
- ✅ Built-in sliding window analysis
- ✅ Automatic consistent rooting
- ✅ Publication-quality recombination detection
- ⚠️ Commercial license required

**SlidingBayes** (Oxford Bioinformatics)
- ✅ Bayesian sliding window framework
- ✅ Natural uncertainty handling
- ✅ Bootstrap support visualization
- ✅ Free academic use

## 📊 Validation and Assessment

### Rooting Consistency Validation

```bash
# Validate rooting across all trees
python validate_rooting_consistency.py results_*/best_rooted_trees.newick
```

**Metrics Assessed:**
- Root depth consistency
- Sister group stability
- Split pattern conservation
- Overall consistency score (target: >0.95)

### Topological Congruence (TreeCmp)

```bash
# Comprehensive tree comparison
treecmp -i tree_list.txt -o comparison_results.txt -d rf ms pd
```

**Metrics Available:**
- Robinson-Foulds distance
- Matching split distance
- Path difference metrics
- Transfer distance

## 🎯 Expected Improvements

Based on literature and tool capabilities:

| Metric                | Current  | Target | Method              |
| --------------------- | -------- | ------ | ------------------- |
| Mean RF Distance      | 0.87     | <0.50  | Constraint trees    |
| Root Consistency      | ~60%     | >95%   | IQ-TREE constraints |
| Topological Stability | Moderate | High   | Professional tools  |

## 📁 File Organization

```
automated-window-sliding/
├── ADVANCED_ROOTING_TOOLS_SURVEY_2024.md      # Comprehensive tool survey
├── advanced_rooting_pipeline.py                # Full-featured pipeline
├── implement_constraint_rooting.py             # IQ-TREE constraint implementation
├── validate_rooting_consistency.py             # Rooting validation tool
├── run_constraint_test.py                      # Quick test script
├── results_constraint_test/                    # Test results
├── results_improved_strategy2/                 # Previous results
└── ingroup_only_alignment.fasta               # Input data
```

## 🔬 Strategy Selection Guide

### For Immediate Improvement (Today)
```bash
python run_constraint_test.py
```
**Best for:** Quick validation of constraint-based approach

### For Research/Development
```bash
python advanced_rooting_pipeline.py --strategy constraint_guided
```
**Best for:** Testing multiple approaches and parameter optimization

### For Publication Quality
- **RDP4/RDP5**: Professional recombination analysis with built-in consistent rooting
- **SlidingBayes**: Bayesian framework with uncertainty quantification
- **IQ-TREE + TreeCmp**: Custom pipeline with comprehensive validation

## 🎯 Success Criteria

### Excellent Results (>95% consistency)
- ✅ Use results for publication
- ✅ Document methodology
- ✅ Consider this the final approach

### Good Results (85-95% consistency)
- ✅ Acceptable for most analyses
- 🔧 Minor optimization may help
- 📊 Compare with professional tools

### Moderate Results (70-85% consistency)
- 🔧 Try non-reversible models
- 🔧 Adjust window parameters
- 🔧 Consider RDP4/RDP5

### Poor Results (<70% consistency)
- 🔧 **Strongly recommended:** Professional tools
- 🔧 Evaluate alignment quality
- 🔧 Consider larger windows or different step sizes

## 📚 Key References

- **IQ-TREE Rooting**: Naser-Khdour et al. (2021) - Non-reversible models for root placement
- **RAxML-NG**: Kozlov et al. (2019) - Constraint trees and topology optimization
- **TreeCmp**: Bogdanowicz et al. (2012) - Comprehensive phylogenetic tree comparison
- **SlidingBayes**: Oxford Bioinformatics - Bayesian sliding window analysis
- **Visual TreeCmp**: Goluch (2020) - Web-based tree comparison platform

## 🆘 Troubleshooting

### IQ-TREE Not Found
```bash
# Install via conda
conda install -c bioconda iqtree

# Or via package manager
sudo apt-get install iqtree  # Ubuntu/Debian
brew install iqtree         # macOS
```

### Memory Issues
- Reduce bootstrap replicates: `-B 100` instead of `-B 1000`
- Use fewer threads: `-T 1` instead of `-T AUTO`
- Process windows in smaller batches

### Poor Consistency Results
1. **Check alignment quality** - Remove problematic sequences
2. **Adjust window parameters** - Try 300-500 bp windows
3. **Evaluate outgroup** - Ensure appropriate evolutionary distance
4. **Consider professional tools** - RDP4/RDP5 for publication quality

## 🎉 Quick Win Strategy

For immediate results on your norovirus dataset:

1. **Run the test** (5 minutes setup):
   ```bash
   python run_constraint_test.py
   ```

2. **Check consistency score** (in output):
   - >0.95: Excellent! Use these results
   - 0.85-0.95: Good! Minor optimization possible
   - <0.85: Try non-reversible models or professional tools

3. **Use improved trees**:
   - Trees saved in `results_constraint_test/all_constraint_trees.newick`
   - Validation report in `results_constraint_test/rooting_validation/`

## 🔮 Future Directions

- **Integration with Nextflow pipeline**: Automate constraint-based rooting
- **Parameter optimization**: Systematic testing of window sizes and models
- **Professional tool evaluation**: RDP4/RDP5 comparison study
- **Publication preparation**: Document methodology and validation results

---

**Status**: Ready for immediate implementation and testing
**Confidence**: High - based on comprehensive 2024-2025 tool survey
**Next Step**: Run `python run_constraint_test.py` to test on your norovirus dataset
