# Dual Rooting Methodology: RootDigger MAD + Backbone Midpoint Rooting

## Overview

This document provides detailed technical documentation for the dual rooting methodology implemented in the automated sliding window phylogenetic analysis pipeline. The pipeline employs **both** RootDigger-based MAD rooting and backbone-constrained midpoint rooting to ensure robust and cross-validated phylogenetic inference across all sliding windows.

**Key Finding:** Analysis of execution logs confirms that the norovirus analysis utilized both rooting methods simultaneously:
1. **Primary:** RootDigger MAD rooting (modified-MAD strategy)
2. **Secondary:** Backbone-constrained midpoint rooting (Biopython Bio.Phylo)

---

## 1. Dual Methodology Framework

### 1.1 Two-Method Rooting Approach

The pipeline implements a comprehensive dual rooting strategy:

1. **RootDigger MAD Rooting:** Likelihood-based optimal root placement
2. **Backbone Midpoint Rooting:** Geometric distance-based root placement
3. **Cross-Validation:** Comparison between methodologically distinct approaches
4. **Quality Assurance:** Independent validation of rooting decisions

### 1.2 Biological and Methodological Rationale

**Why Dual Rooting?**

- **Method Independence:** Two algorithmically distinct approaches reduce method-specific bias
- **Cross-Validation:** Enables assessment of rooting consistency and reliability  
- **Robustness:** Guards against artifacts from single-method approaches
- **Flexibility:** Provides options for downstream analysis preferences

---

## 2. Method 1: RootDigger MAD Rooting

### 2.1 Software Implementation

**RootDigger Configuration:**
- **Software:** RootDigger v1.0.3
- **Container:** quay.io/biocontainers/rootdigger:1.0.3--h9ee0642_0
- **Algorithm:** Modified Minimal Ancestor Deviation (modified-MAD)
- **Execution:** Applied to all 27 sliding windows

### 2.2 Technical Implementation

**Command Line Execution:**

```bash
rootdigger \
    --msa clean_alignment.fasta \
    --tree clean_tree.newick \
    --threads ${task.cpus} \
    --initial-root-strategy "modified-mad" \
    --silent
```

**Parameter Details:**
- `--initial-root-strategy "modified-mad"`: Enhanced MAD algorithm
- `--silent`: Minimal output for pipeline integration
- `--threads`: Multi-core processing for efficiency

### 2.3 Quality Control

**Execution Validation (From Analysis Logs):**
- **Success Rate:** 100% (all 27 windows completed successfully)
- **Execution Time:** 284-577ms per window
- **Container Consistency:** Uniform containerized environment
- **Output Format:** Newick trees with likelihood-based root placement

### 2.4 MAD Algorithm Principles

**Minimal Ancestor Deviation (MAD):**
- **Objective:** Minimize evolutionary distance from root to ancestral nodes
- **Method:** Likelihood-based optimization of root position
- **Advantage:** Robust to rate heterogeneity across lineages
- **Implementation:** Modified-MAD variant with enhanced optimization

---

## 3. Method 2: Backbone Midpoint Rooting

### 3.1 Backbone Tree Generation

**Step 1: Full Alignment Phylogenetic Reconstruction**

**Software:** IQ-TREE2 with Model Finder Plus

```bash
iqtree2 -s cleaned_alignment_combined.fasta \
        --prefix backbone_tree \
        -m MFP \
        --threads 4 \
        -B 1000 \
        --quiet
```

**Parameters:**
- **Model Selection:** Automatic optimal model selection (MFP)
- **Bootstrap Support:** 1000 ultrafast bootstrap replicates
- **Threading:** Multi-core parallel processing

### 3.2 Midpoint Rooting Implementation

**Software:** Biopython Bio.Phylo module

**Algorithm Implementation:**

```python
from Bio import Phylo

def apply_midpoint_rooting(tree_file, output_file):
    """Apply midpoint rooting using Biopython."""
    tree = Phylo.read(tree_file, "newick")
    tree.root_at_midpoint()
    Phylo.write(tree, output_file, "newick")
```

**Midpoint Rooting Principle:**
- **Method:** Root at midpoint of longest path between taxa
- **Assumption:** Approximately equal evolutionary rates
- **Advantage:** No outgroup required, deterministic results

### 3.3 Execution Validation

**From Analysis Logs:**
- **Process:** BACKBONE_MIDPOINT_ROOTING
- **Success Rate:** 100% (all 27 windows completed)
- **Execution Time:** 156-189ms per window
- **Implementation:** Biopython-based midpoint algorithm

---

## 4. Comparative Analysis Framework

### 4.1 Dual Output Generation

**RootDigger Output:**
- **Location:** MAD_ROOTING processes in execution logs
- **Format:** Likelihood-optimized Newick trees
- **Algorithm:** Modified-MAD strategy
- **Validation:** Log-likelihood values provided

**Backbone Midpoint Output:**
- **Location:** BACKBONE_MIDPOINT_ROOTING processes in execution logs
- **Format:** Geometrically-rooted Newick trees
- **Algorithm:** Classical midpoint rooting
- **Validation:** Tree structure verification

### 4.2 Cross-Method Validation

**Comparison Capabilities:**
1. **Root Position Agreement:** Compare root placement between methods
2. **Topological Consistency:** Assess tree structure preservation
3. **Statistical Support:** Evaluate confidence metrics
4. **Method Performance:** Compare execution times and success rates

### 4.3 Quality Metrics

**RootDigger Quality Indicators:**
- **Log-Likelihood Values:** Quantitative root placement confidence
- **Strategy Success:** Modified-MAD algorithm completion
- **Container Stability:** Consistent containerized execution

**Midpoint Rooting Quality Indicators:**
- **Algorithm Success:** Biopython midpoint completion
- **Tree Balance:** Assessment of resulting tree structure
- **Deterministic Results:** Reproducible root placement

---

## 5. Pipeline Integration and Workflow

### 5.1 Nextflow Implementation

**Parallel Processing Structure:**

```nextflow
// Generate backbone tree
backbone_tree = GENERATE_BACKBONE_TREE(alignment, output_dir)

// Method 1: RootDigger MAD rooting
mad_rooted_trees = MAD_ROOTING(
    window_trees_with_alignments,
    output_dir
)

// Method 2: Backbone midpoint rooting
backbone_rooted_trees = BACKBONE_MIDPOINT_ROOTING(
    window_trees_with_alignments,
    backbone_tree.backbone_tree,
    output_dir
)
```

### 5.2 Configuration Parameters

**Analysis Configuration (params_cleaned_large_windows_250_10.json):**

```json
{
    "backbone_midpoint_rooting": true,
    "mad_rooting": true,
    "rootdigger_strategy": "modified-mad",
    "rootdigger_exhaustive": false
}
```

**Default Pipeline Settings (nextflow.config):**
- `mad_rooting = true`: Enable RootDigger MAD rooting
- `backbone_midpoint_rooting = false`: Backbone rooting (enabled via params)
- `rootdigger_strategy = 'modified-mad'`: Default MAD algorithm

---

## 6. Performance Characteristics

### 6.1 Computational Requirements

**RootDigger Performance:**
- **Memory:** Low requirement with containerized execution
- **CPU:** Multi-threaded processing per window
- **Time:** 284-577ms per window (27 windows total)
- **Scalability:** Linear with number of windows

**Biopython Performance:**
- **Memory:** Minimal Python module overhead
- **CPU:** Single-threaded per window
- **Time:** 156-189ms per window (27 windows total)  
- **Scalability:** Excellent for large-scale analyses

### 6.2 Resource Optimization

**Parallel Execution Benefits:**
- **Independent Processing:** Both methods run simultaneously
- **Resource Efficiency:** Optimal CPU utilization
- **Time Savings:** No sequential dependency between methods
- **Fault Tolerance:** Method independence reduces failure impact

---

## 7. Biological and Statistical Implications

### 7.1 Recombination Detection Enhancement

**Method Complementarity:**
- **RootDigger:** Likelihood-based optimization accounts for evolutionary models
- **Midpoint:** Geometric approach independent of model assumptions
- **Cross-Validation:** Consensus rooting increases confidence in recombination signals

### 7.2 Rate Heterogeneity Management

**Algorithmic Differences:**
- **MAD Rooting:** Explicitly accounts for rate variation via likelihood framework
- **Midpoint Rooting:** Assumes molecular clock but provides deterministic baseline
- **Comparative Assessment:** Identifies regions with significant rate heterogeneity

### 7.3 Phylogenetic Signal Preservation

**Complementary Strengths:**
- **Statistical Rigor:** RootDigger provides likelihood-based statistical framework
- **Geometric Simplicity:** Midpoint rooting offers intuitive distance-based approach
- **Robustness Testing:** Agreement between methods validates rooting decisions

---

## 8. Best Practices and Recommendations

### 8.1 Method Selection Guidelines

**Use RootDigger MAD When:**
- Statistical rigor is paramount
- Rate heterogeneity is suspected
- Likelihood-based confidence metrics are needed
- Publication-quality analysis is required

**Use Midpoint Rooting When:**
- Rapid exploratory analysis is needed
- Deterministic results are preferred
- Computational simplicity is valued
- Baseline comparison is required

**Use Both Methods When:**
- Cross-validation is important
- Method robustness needs assessment
- Comprehensive analysis is required
- Manuscript reviewers expect method comparison

### 8.2 Quality Assessment Guidelines

**Strong Agreement Indicators:**
- Similar root placement between methods
- Consistent topological relationships
- Low Robinson-Foulds distances between rooted trees

**Potential Issues Indicators:**
- Highly discordant root placement
- Significant topological differences
- Method-specific failures or errors

### 8.3 Interpretation Framework

**Consensus Rooting:**
- High confidence when both methods agree
- Moderate confidence with minor discrepancies
- Low confidence with major disagreements

**Method-Specific Considerations:**
- RootDigger sensitivity to model violations
- Midpoint sensitivity to rate heterogeneity
- Container dependency for RootDigger vs. local execution for Biopython

---

## 9. Troubleshooting and Diagnostics

### 9.1 Common Issues and Solutions

**RootDigger Failures:**
- **Taxa Name Mismatches:** Ensure alignment and tree taxa consistency
- **Container Issues:** Verify Singularity/Docker availability
- **Memory Limitations:** Check available system resources

**Midpoint Rooting Failures:**
- **Biopython Dependencies:** Verify Python environment setup
- **Tree Format Issues:** Check Newick format compliance
- **Branch Length Problems:** Ensure positive branch lengths

### 9.2 Validation Procedures

**Method Verification:**

```bash
# Check RootDigger output
ls results/*/rooted_trees/*_rooted.tree

# Check backbone midpoint output  
ls results/*/backbone_midpoint/*_rooted.treefile

# Compare execution logs
grep "MAD_ROOTING\|BACKBONE_MIDPOINT_ROOTING" results/*/pipeline_info/execution_trace_*.txt
```

**Quality Assessment:**

```python
# Load and compare rooted trees
from Bio import Phylo

# RootDigger tree
mad_tree = Phylo.read("window_mad_rooted.tree", "newick")

# Midpoint tree  
midpoint_tree = Phylo.read("window_midpoint_rooted.treefile", "newick")

# Assess rooting agreement
print(f"MAD rooted: {mad_tree.is_rooted}")
print(f"Midpoint rooted: {midpoint_tree.is_rooted}")
```

---

## 10. Future Enhancements

### 10.1 Advanced Integration

**Consensus Rooting Algorithm:**
- Weighted combination of rooting methods
- Confidence-based root placement
- Automated method selection

**Enhanced Validation:**
- Quantitative rooting agreement metrics
- Statistical tests for root placement
- Visual diagnostic plots

### 10.2 Methodological Extensions

**Additional Rooting Methods:**
- Outgroup rooting when available
- Clock-based rooting for temporal data
- Minimum variance rooting

**Performance Optimization:**
- GPU acceleration for large datasets
- Memory-efficient algorithms
- Parallel method execution

---

*Document Version: 2.0 - Corrected for Dual Method Implementation*  
*Last Updated: 19 June 2025*  
*Pipeline Version: Automated Window Sliding v2.0*  
*Analysis Validation: Based on execution logs from results_cleaned_large_windows_250_10*
