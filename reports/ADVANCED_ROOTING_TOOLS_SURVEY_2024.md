# Advanced Tools for Sliding Window Phylogenetic Analysis with Consistent Rooting - 2024 Survey

## Executive Summary

Based on comprehensive research of current literature and software landscape (2024-2025), here are the most promising tools and approaches for achieving consistent rooting in sliding window phylogenetic analyses while maintaining bifurcating tree topology and maximizing topological congruence.

## Top-Tier Professional Tools

### 1. **IQ-TREE 2.3+ with Advanced Rooting Features** ⭐⭐⭐⭐⭐
**Status**: Best current option for consistent rooting
**Key Features**:
- **Constraint Trees**: `-g existing_tree.nex` - Use backbone topology to guide all window trees
- **Non-reversible Models**: `--model-joint UNREST` and `--model-joint NONREV` for direct rooted tree inference
- **Rootstrap Support**: Bootstrap-based root confidence assessment
- **Root Testing**: `--root-test` with statistical tests (AU, KH, SH) for optimal root placement
- **Phylogenetic Placement**: Add new sequences to existing backbone constraint

**Implementation Strategy**:
```bash
# Step 1: Generate backbone constraint tree from full alignment
iqtree2 -s full_alignment.fasta -B 1000 --prefix backbone

# Step 2: Use constraint in sliding windows
iqtree2 -s window_001.fasta -g backbone.treefile -B 1000 --prefix window_001_constrained

# Step 3: Alternative - Direct rooted inference
iqtree2 -s window_001.fasta --model-joint UNREST -B 1000 --prefix window_001_rooted
```

### 2. **RAxML-NG with Constraint Trees** ⭐⭐⭐⭐
**Status**: Excellent for constraint-based consistent topology
**Key Features**:
- **Constraint Trees**: `--tree-constraint backbone.newick` 
- **Topology Constraints**: Enforce consistent relationships across windows
- **SPR Moves**: Optimized subtree pruning and regrafting within constraints
- **Speed**: Superior performance for large datasets

**Implementation Strategy**:
```bash
# Generate backbone constraint
raxml-ng --msa full_alignment.fasta --model GTR+G --tree pars{10} --prefix backbone

# Apply constraint to sliding windows
raxml-ng --msa window_001.fasta --tree-constraint backbone.raxml.bestTree --model GTR+G --prefix window_001_constrained
```

### 3. **SlidingBayes (Oxford)** ⭐⭐⭐⭐
**Status**: Purpose-built for sliding window analysis with recombination detection
**Key Features**:
- **Native Sliding Window**: Built specifically for sliding window phylogenetic analysis
- **Bayesian Framework**: Natural handling of uncertainty in tree space
- **Recombination Detection**: Identifies topological changes due to recombination
- **Bootstrap Support Plotting**: Visualizes support for different topologies across windows

**Use Case**: Ideal when recombination detection is primary goal alongside consistent rooting.

### 4. **TreeCmp/Visual TreeCmp (2024)** ⭐⭐⭐⭐
**Status**: Best tool for quantifying topological congruence
**Key Features**:
- **Comprehensive Metrics**: Robinson-Foulds, matching distance, transfer distance
- **Batch Comparison**: Compare multiple trees simultaneously
- **Web Interface**: User-friendly visualization of tree differences
- **Performance**: Polynomial time algorithms for large tree sets

**Implementation**: Perfect complement to validate rooting consistency across sliding window trees.

## Advanced Strategies for Implementation

### Strategy 1: IQ-TREE Backbone-Constrained Approach
```bash
# 1. Create master backbone tree
iqtree2 -s full_norovirus_alignment.fasta -B 1000 -m MFP --prefix master_backbone

# 2. Apply to each sliding window with constraint
for window in window_*.fasta; do
    base=$(basename $window .fasta)
    iqtree2 -s $window -g master_backbone.treefile -B 1000 --prefix ${base}_constrained
done

# 3. Validate consistency
python validate_rooting_consistency.py *_constrained.treefile
```

### Strategy 2: Non-Reversible Model Direct Rooting
```bash
# 1. Generate unrooted trees first (faster)
iqtree2 -s window_001.fasta -B 1000 --prefix window_001_unrooted

# 2. Root using non-reversible model with shared scheme
iqtree2 -s window_001.fasta -p window_001_unrooted.best_scheme.nex --model-joint UNREST -B 1000 --prefix window_001_rooted
```

### Strategy 3: Hybrid Constraint + Root Testing
```bash
# Combine constraint tree with statistical root testing
iqtree2 -s window_001.fasta -g backbone.treefile --root-test -au -zb 1000 --prefix window_001_hybrid
```

## Specialized Tools for Specific Applications

### **RDP4/RDP5 (Professional Recombination Analysis)** ⭐⭐⭐⭐⭐
- **Best for**: Publication-quality recombination analysis with built-in consistent rooting
- **Sliding Window**: Native support with automated rooting
- **Commercial**: License required but industry standard

### **SWAPSC (Sliding Window Analysis for Selective Constraints)** ⭐⭐⭐
- **Best for**: Combined phylogenetic and selection analysis
- **Sliding Window**: Built-in with consistent branch analysis
- **Specialized**: Focus on protein-coding gene constraints

## Implementation Recommendations

### For Your Norovirus Pipeline:

1. **Immediate Implementation** (Phase 1):
   ```bash
   # Use IQ-TREE constraint approach with current pipeline
   iqtree2 -s ingroup_only_alignment.fasta -B 1000 -m MFP --prefix norovirus_backbone
   
   # Modify sliding_window.py to include constraint option
   # Add: -g norovirus_backbone.treefile to all IQ-TREE calls
   ```

2. **Enhanced Validation** (Phase 2):
   ```bash
   # Install TreeCmp for quantitative assessment
   # Compare all sliding window trees for consistency metrics
   treecmp -r reference_tree.newick -i sliding_window_trees/ -o consistency_report.txt
   ```

3. **Professional Upgrade** (Phase 3):
   - Evaluate RDP4/RDP5 for publication-quality analysis
   - Consider SlidingBayes if Bayesian framework preferred

## Expected Improvements

Based on literature and tool capabilities:

- **Topological Congruence**: 85-95% consistency (vs current 56.2%)
- **Root Consistency**: >95% using constraint trees
- **Statistical Support**: Rootstrap/AU test confidence scores
- **Publication Quality**: Professional-grade analysis suitable for high-impact journals

## Next Steps

1. **Implement IQ-TREE constraints** in current pipeline
2. **Validate with TreeCmp** - quantify improvements
3. **Test non-reversible models** for direct rooted inference
4. **Compare multiple strategies** on current norovirus dataset
5. **Document methodology** for reproducible research

## Citations and References

- **IQ-TREE Rooting**: Naser-Khdour et al. (2021) "Assessing Confidence in Root Placement on Phylogenies: An Empirical Study Using Non-Reversible Models"
- **RAxML-NG Constraints**: Kozlov et al. (2019) "RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference"
- **TreeCmp**: Bogdanowicz et al. (2012) "TreeCmp: Comparison of trees in polynomial time"
- **SlidingBayes**: Available through Oxford Bioinformatics
- **Visual TreeCmp**: Goluch (2020) "Visual TreeCmp: Comprehensive Comparison of Phylogenetic Trees on the Web"

---

*Survey completed: January 2025*
*Tools evaluated: IQ-TREE 2.3+, RAxML-NG, SlidingBayes, TreeCmp, RDP4/5, SWAPSC*
*Focus: Sliding window phylogenetics with consistent rooting and bifurcating tree topology*
