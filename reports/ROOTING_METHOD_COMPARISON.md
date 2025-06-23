# Sliding Window Phylogenetic Analysis - Rooting Method Comparison

**Generated:** 2025-06-18T16:23:55.669145

## Overview

This document compares two optimized approaches for rooting phylogenetic trees in sliding window analysis:

1. **RootDigger Exhaustive Rooting** (Fully optimized)
2. **Backbone-Constrained Midpoint Rooting** (Fast with backbone guidance)

## Analysis Parameters (Both Methods)

- **Sequence Data:** 38 Norovirus sequences (refined alignment)
- **Window Size:** 200 bp
- **Step Size:** 20 bp
- **Total Windows:** 25
- **Backbone Constraint:** IQ-TREE backbone tree used for guidance
- **Taxa Name Cleaning:** Applied for compatibility

## Method 1: RootDigger Exhaustive Rooting

### Approach
- **Tree Reconstruction:** IQ-TREE with backbone constraint
- **Rooting Method:** RootDigger exhaustive search
- **Optimization Level:** Maximum (exhaustive root testing)
- **Speed:** Slower due to exhaustive optimization

### Results
- **Success Rate:** 100% (25/25 windows)
- **Root Quality:** Optimal root placement based on likelihood
- **Binary Roots:** 100% validated
- **Output Files:** 
  - `ALL_ROOTED_TREES_OPTIMIZED.tre` (37.6 KB)
  - `ALL_ROOTED_TREES_SIMPLE.tre` (34.5 KB)

### Advantages
- ✅ **Optimal root placement** through exhaustive search
- ✅ **Maximum likelihood rooting** for each window
- ✅ **Publication-quality** results
- ✅ **Best for detailed phylogenetic analysis**

### Disadvantages
- ⏱️ **Slower execution** (exhaustive search per tree)
- 🧮 **Computationally intensive**

## Method 2: Backbone-Constrained Midpoint Rooting

### Approach
- **Tree Reconstruction:** IQ-TREE with backbone constraint
- **Rooting Method:** Midpoint rooting (Biopython)
- **Optimization Level:** Fast with backbone guidance
- **Speed:** Much faster execution

### Results
- **Success Rate:** 100% (25/25 windows)
- **Root Quality:** Balanced rooting based on branch lengths
- **Binary Roots:** 100% validated
- **Output Files:**
  - `ALL_BACKBONE_MIDPOINT_TREES.tre` (annotated)
  - `ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre` (simple)

### Advantages
- ⚡ **Fast execution** (quick midpoint calculation)
- 🎯 **Backbone constraint optimization** maintained
- ✅ **Reliable rooting** for comparative analysis
- 📊 **Good for large-scale studies**

### Disadvantages
- 🔬 **Less optimal** than exhaustive search
- 📏 **Depends on branch length accuracy**

## Performance Comparison

| Metric | RootDigger Exhaustive | Backbone Midpoint |
|--------|----------------------|-------------------|
| **Execution Time** | ~10-15 min per tree | ~1-2 min per tree |
| **Root Optimality** | Maximum likelihood | Branch length based |
| **Backbone Constraint** | ✅ Applied | ✅ Applied |
| **Success Rate** | 100% | 100% |
| **Binary Roots** | 100% | 100% |
| **Best For** | Final analysis | Preliminary/large studies |

## File Structure Comparison

### RootDigger Exhaustive Results
```
results_optimized_200_20/
├── backbone_analysis/          # Backbone constraint tree
├── constraint_trees/           # Backbone-constrained trees
├── rooted_trees_exhaustive/    # RootDigger optimized trees
├── sliding_windows/            # Window alignments
└── ANALYSIS_SUMMARY.md         # Comprehensive report
```

### Backbone Midpoint Results
```
results_backbone_midpoint_200_20/
├── backbone_tree.treefile      # Reference backbone
├── constraint_trees_midpoint/  # Backbone-constrained trees
└── rooted_trees_midpoint/      # Midpoint rooted trees
```

## Recommendations

### Use RootDigger Exhaustive When:
- 🎯 **Publication-quality analysis** is required
- 🔬 **Detailed phylogenetic relationships** are crucial
- ⏱️ **Time is not a constraint**
- 📊 **Small to medium datasets** (≤100 trees)

### Use Backbone Midpoint When:
- ⚡ **Speed is important**
- 📈 **Large-scale comparative studies**
- 🔄 **Preliminary analysis** before detailed rooting
- 💻 **Computational resources are limited**

## Tree Quality Assessment

Both methods produce high-quality results with:
- ✅ **Proper binary root structures**
- ✅ **Backbone constraint optimization**
- ✅ **Taxa name consistency**
- ✅ **Valid Newick format trees**
- ✅ **100% success rates**

## Combined Analysis Benefits

The backbone constraint optimization in both methods provides:

1. **Phylogenetic Consistency:** All window trees guided by backbone topology
2. **Reduced Tree Space:** Constraint reduces search complexity
3. **Improved Accuracy:** Backbone prevents spurious relationships
4. **Computational Efficiency:** Even exhaustive search is faster with constraints

## Conclusion

Both approaches successfully generate high-quality rooted phylogenetic trees:

- **RootDigger Exhaustive:** Best for final, publication-quality analysis
- **Backbone Midpoint:** Excellent for fast, reliable comparative studies

The choice depends on your specific needs for optimality vs. speed, with both methods benefiting from backbone constraint optimization.

---

## Output Files Available

### RootDigger Exhaustive
- `ALL_ROOTED_TREES_OPTIMIZED.tre` - Fully annotated trees
- `ALL_ROOTED_TREES_SIMPLE.tre` - Simple Newick format
- `TREE_SUMMARY_TABLE.txt` - Window position index

### Backbone Midpoint  
- `ALL_BACKBONE_MIDPOINT_TREES.tre` - Fully annotated trees
- `ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre` - Simple Newick format
- `BACKBONE_MIDPOINT_TREE_SUMMARY.txt` - Window position index

*Both sets ready for phylogenetic analysis software!*
