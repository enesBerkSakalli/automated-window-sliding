# RootDigger Optimization Strategies

This document outlines the best practices for using RootDigger in the Automated Window Sliding pipeline based on current research and our implementation.

## Overview

RootDigger is a phylogenetic rooting tool that uses non-reversible Markov models to find the most likely root location on a tree. It's particularly useful for our sliding window analysis because it:

- Works without requiring outgroups
- Provides consistent rooting across window segments
- Is faster than IQ-TREE for most datasets
- Can quantify uncertainty in root placement

## Current Implementation

### Default Strategy: `modified-mad`
- **Best for**: Most general use cases, balanced speed vs accuracy
- **Advantages**: Improved MAD approach, faster than exhaustive search
- **Use when**: Running large-scale analyses or when computational time is a constraint

### Alternative Strategies

#### 1. Exhaustive Mode (`--exhaustive`)
- **Best for**: High-accuracy requirements, smaller datasets
- **Research finding**: "Exhaustive mode is more successful at identifying the correct root location"
- **Trade-off**: Significantly slower but more thorough search
- **Use when**: Root placement accuracy is critical, dataset size is manageable

#### 2. Traditional MAD (`mad-like`)
- **Best for**: Comparison with classic MAD implementations
- **Use when**: Benchmarking against traditional methods

#### 3. Midpoint Rooting (`midpoint`)
- **Best for**: Quick approximation, molecular clock assumption holds
- **Use when**: Fast preliminary analysis needed

## Configuration Options

### Basic Configuration
```bash
# Use default modified-MAD strategy (recommended)
nextflow run main.nf --input data/alignment.fasta

# Use exhaustive mode for higher accuracy (slower)
nextflow run main.nf --input data/alignment.fasta --rootdigger_exhaustive true

# Use traditional MAD strategy (deprecated - use random instead)
nextflow run main.nf --input data/alignment.fasta --rootdigger_strategy random
```

### Advanced Configuration
```bash
# For high-accuracy analysis with exhaustive search
nextflow run main.nf \
    --input data/alignment.fasta \
    --rootdigger_strategy modified-mad \
    --rootdigger_exhaustive true \
    --window_size 300 \
    --step_size 50

# For fast analysis with midpoint rooting
nextflow run main.nf \
    --input data/alignment.fasta \
    --rootdigger_strategy midpoint \
    --window_size 1000 \
    --step_size 200
```

## Performance Considerations

### Dataset Size vs Strategy
- **Small datasets (< 50 taxa)**: Consider exhaustive mode for best accuracy
- **Medium datasets (50-200 taxa)**: Modified-MAD provides good balance
- **Large datasets (> 200 taxa)**: Modified-MAD or random for efficiency

### Window Analysis Specific
- **Sliding windows benefit from consistent rooting**: MAD-based strategies preferred
- **Multiple windows**: Speed becomes important, avoid exhaustive unless necessary
- **Viral phylogenetics**: Modified-MAD works well for fast-evolving sequences

## Research-Based Recommendations

### When to Use Exhaustive Mode
Based on research findings that "exhaustive mode is more successful at identifying the correct root location":

1. **High-stakes analyses** where root placement accuracy is critical
2. **Datasets with conflicting phylogenetic signal**
3. **Benchmarking or validation studies**
4. **Small to medium datasets** where runtime is manageable

### When to Use Standard Mode
For routine sliding window analyses:

1. **Large-scale genomic studies**
2. **Exploratory analyses**
3. **Production pipelines** where speed matters
4. **Comparative studies** across many datasets

## Quality Assessment

The pipeline logs rooting confidence information when available:
- Check `*_rooting.log` files for strategy used
- Look for log-likelihood values indicating root support
- Compare results across different strategies if needed

## Troubleshooting

### Common Issues
1. **Taxa name mismatches**: Pipeline automatically cleans headers
2. **Failed rooting**: Falls back gracefully with error logging
3. **Memory issues**: Reduce to non-exhaustive mode for large datasets

### Validation
- Compare rooted trees across strategies
- Check for consistent root placement across windows
- Verify biological plausibility of root locations

## Future Improvements

Potential enhancements being considered:
- Automatic strategy selection based on dataset characteristics
- Parallel exhaustive search optimization
- Integration with other rooting methods for comparison
- Automated quality metrics for root placement confidence

## References

- Bettisworth, B., & Stamatakis, A. (2021). RootDigger: a root placement program for phylogenetic trees. BMC Bioinformatics, 22(1), 1-13.
- Tria, F. D., Landan, G., & Dagan, T. (2017). Phylogenetic rooting using minimal ancestor deviation. Nature Ecology & Evolution, 1(7), 1-7.
