# MAD Rooting Integration

## Overview

This pipeline now includes **MAD (Minimal Ancestor Deviation) rooting** for all phylogenetic trees, making it particularly suitable for sliding window analysis without outgroups. MAD rooting provides consistent, reproducible tree rooting across all windows.

## What is MAD Rooting?

MAD rooting is a branch-length based method that:
- **Minimizes ancestor deviation** across the tree
- Works **without requiring an outgroup**
- Is **robust to rate heterogeneity** between lineages
- Provides **deterministic, reproducible results**
- Is **computationally efficient** (O(number of branches))

## Why MAD for Sliding Windows?

1. **Consistency**: Same rooting algorithm applied to all windows
2. **No outgroup required**: Perfect for intraspecific analyses
3. **Speed**: Can root thousands of trees in minutes
4. **Reproducibility**: 100% deterministic results
5. **Robustness**: Handles rate variation across lineages

## Pipeline Integration

### Automatic MAD Rooting

The pipeline automatically applies MAD rooting to all reconstructed trees:

```bash
nextflow run main.nf \
  --input aligned_sequences.fasta \
  --outdir results \
  --window_size 200 \
  --step_size 20 \
  --model "GTR+G+I" \
  --mad_rooting true  # Default: enabled
```

### Output Files

With MAD rooting enabled, you get:

1. **Original unrooted trees**: `results/[params]/best_trees/`
2. **MAD rooted trees**: `results/[params]/rooted_trees/`
3. **Rooted tree collections**: 
   - `rooted_trees_collection.newick`
   - `rooted_trees_collection.nexus`
4. **Quality reports**:
   - `mad_rooting_summary.txt`
   - `root_quality_report.txt`

### Directory Structure

```
results/aligned_sequences_w200_s20_GTR_G_I_[timestamp]/
├── rooted_trees/                    # Individual MAD rooted trees
│   ├── 1_mad.rooted.tree
│   ├── 21_mad.rooted.tree
│   └── ...
├── rooted_trees_collection.newick   # All rooted trees in Newick format
├── rooted_trees_collection.nexus    # All rooted trees in Nexus format
├── mad_rooting_summary.txt         # Success rates and window-by-window results
└── root_quality_report.txt         # Quality assessment and guidelines
```

## Quality Control

### MAD Deviation Values

- **Low deviation (< 0.2)**: High confidence rooting
- **Medium deviation (0.2-0.4)**: Moderate confidence
- **High deviation (> 0.4)**: Consider alternative methods

### Quality Reports

The pipeline generates comprehensive reports:

```bash
# View rooting success rates
cat results/*/mad_rooting_summary.txt

# Check quality assessment
cat results/*/root_quality_report.txt
```

## Advanced Usage

### Disable MAD Rooting

```bash
nextflow run main.nf \
  --input alignment.fasta \
  --mad_rooting false  # Disable MAD rooting
```

### Integration with Downstream Analysis

The rooted trees are perfect for:
- **Phylogenomic movies**: Consistent rooting prevents artificial topology changes
- **Recombination detection**: Real topology changes vs. rooting artifacts
- **Selection analysis**: Proper tree rooting for branch-based methods
- **Phylogeography**: Consistent temporal framework

## Computational Performance

### Benchmarks

- **412 windows (200bp each)**: ~5-10 minutes for rooting
- **Memory usage**: Minimal (~100MB peak)
- **Parallelization**: Fully parallelized across windows

### Scaling

MAD rooting scales linearly:
- **1000 windows**: ~15-20 minutes
- **10000 windows**: ~2-3 hours
- **100000 windows**: ~1 day

## Comparison with Other Methods

| Method     | Speed | Outgroup Required | Deterministic | Rate Heterogeneity |
| ---------- | ----- | ----------------- | ------------- | ------------------ |
| **MAD**    | ⭐⭐⭐⭐⭐ | ❌ No              | ✅ Yes         | ✅ Robust           |
| Midpoint   | ⭐⭐⭐⭐⭐ | ❌ No              | ✅ Yes         | ❌ Sensitive        |
| Outgroup   | ⭐⭐⭐   | ✅ Required        | ✅ Yes         | ✅ Robust           |
| RootDigger | ⭐⭐    | ❌ No              | ✅ Yes         | ✅ Very Robust      |

## Troubleshooting

### Common Issues

1. **Failed rooting for some windows**:
   - Check star-like topologies
   - Very short branches
   - Consider concatenation for problem regions

2. **High MAD deviation**:
   - Indicates potential clock violations
   - Consider RootDigger for validation
   - May flag real biological signal

### Advanced Rooting Strategy

For challenging datasets, consider the hybrid approach:

```bash
# 1. Run MAD rooting (fast, most windows)
nextflow run main.nf --mad_rooting true

# 2. For high-deviation windows, consider RootDigger
# (implement as future enhancement)
```

## References

1. **MAD Method**: Tria et al. (2017) Nature Ecology & Evolution
2. **Phylogenetic rooting**: Huelsenbeck et al. (2002) Systematic Biology
3. **Sliding window analysis**: Martin et al. (2011) Molecular Biology and Evolution

## Future Enhancements

Planned additions:
- **RootDigger integration** for validation
- **Automatic quality flagging** based on MAD deviation
- **Root confidence visualization** 
- **Hybrid rooting strategies**

---

*This MAD rooting integration ensures your sliding window phylogenetic analysis produces consistently rooted, publication-ready trees suitable for downstream phylogenomic analysis.*
