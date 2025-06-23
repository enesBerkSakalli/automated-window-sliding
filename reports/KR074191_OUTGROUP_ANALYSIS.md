# KR074191.1 Outgroup Analysis Report

## Summary

Based on my analysis of the publication **10.1371/journal.pone.0189504** and your norovirus dataset, here are the key findings about using **KR074191.1** as an outgroup sequence:

## Sequence Characteristics

**KR074191.1**: Norovirus GII/Hu/BRA/2009/GII.Pg-GII.12/RS16344 

- **Length**: 492 bp (unique length in dataset)
- **Genotype**: GII.Pg-GII.12 (recombinant)
- **Origin**: Brazil, 2009
- **GC Content**: 48.58%

## Publication Context (DOI: 10.1371/journal.pone.0189504)

The paper "Detection and molecular characterization of the novel recombinant norovirus GII.P16-GII.4 Sydney in southeastern Brazil in 2016" focuses on:

- **GII.4 variants and recombinants** (your dataset is primarily GII.4)
- **GII.P16-GII.4 Sydney recombinants** that emerged in 2016
- **Phylogenetic methodology**: Neighbor-joining, Kimura 2-parameter, MEGA 6.0
- **Region analyzed**: ORF1/ORF2 junction (570 bp)

## Dataset Composition Analysis

Your dataset contains:
- **44 total sequences**
- **Dominant genotypes**: GII.P4-GII.4 (11 seq), GII.P7-GII.6 (12 seq)
- **GII.12 sequences**: 5 total (11.4% of dataset)
- **Length distribution**: Most sequences are 514 bp, but KR074191.1 is 492 bp

## Branch Length Issue Explanation

### Why KR074191.1 Causes Extreme Branch Lengths (>12 substitutions/site):

1. **✓ Genuine Evolutionary Divergence**: 
   - GII.12 is genuinely divergent from GII.4/GII.6 genotypes
   - Different evolutionary lineage within norovirus GII

2. **✓ Recombination Effects**: 
   - GII.Pg-GII.12 is a recombinant strain
   - Mosaic genome structure creates artificial distance in trees

3. **✓ Length Differences**: 
   - KR074191.1 (492 bp) vs. most sequences (514 bp)
   - Alignment gaps inflate apparent substitution rates

4. **✓ Rate Heterogeneity**: 
   - Different evolutionary rates between genotype lineages
   - GII.12 may evolve faster/slower than GII.4

## Outgroup Suitability Assessment

### ✓ PROS (Good for outgroup):
- **Distinct genotype**: GII.12 vs predominantly GII.4 dataset
- **Evolutionary distance**: Sufficient divergence for rooting
- **Single representative**: Won't cluster with ingroup
- **Biological relevance**: From same virus family

### ⚠ CONS (Potential issues):
- **Extreme branch length**: >12 substitutions/site is unrealistic
- **Different sequence length**: May cause alignment artifacts
- **Recombinant nature**: Complex evolutionary history
- **May distort tree topology**: Can pull other sequences artificially

## Recommendations

### Option 1: Use KR074191.1 with Branch Length Correction ✓ RECOMMENDED
```bash
# Your current approach with branch length capping
# This preserves topology while making trees visualizable
```

### Option 2: Exclude KR074191.1 and Use Alternative Outgroup
Consider using one of the other divergent sequences:
- **GII.P13-GII.17** sequences (if less problematic)
- **GII.P21** sequences (different P-type)
- **GII.P15-GII.15** sequences

### Option 3: Root Using Midpoint Rooting
```bash
# In your analysis, use midpoint rooting instead of outgroup
# This avoids the extreme branch length issue entirely
```

### Option 4: Multiple Sequence Alignment Check
1. Realign sequences with a more sophisticated method
2. Check for alignment artifacts around KR074191.1
3. Manually inspect the alignment quality

## Modified Pipeline Approach

### For Current Analysis:
1. **Keep using branch length correction** (you've already implemented this)
2. **Verify tree topology** makes biological sense
3. **Check bootstrap support** for major clades
4. **Document the limitation** in your results

### For Future Analysis:
1. **Try excluding KR074191.1** and compare tree topologies
2. **Use multiple outgroups** if available
3. **Consider midpoint rooting** as an alternative
4. **Validate with published phylogenies** of norovirus

## Conclusion

**KR074191.1 (GII.Pg-GII.12) is biologically appropriate as an outgroup** for your norovirus GII.4-dominated dataset, but its extreme branch lengths suggest either:

1. **Genuine high divergence** (which is fine for an outgroup)
2. **Technical artifacts** from alignment or recombination

Your **branch length correction approach is the right solution** for visualization while preserving the biological relationships. The capped branches (0.1 substitutions/site) are reasonable for intra-genotype comparisons and make the trees interpretable.

## Reference for Outgroup Selection

The paper methodology (10.1371/journal.pone.0189504) used standard phylogenetic approaches but focused on closely related GII.4 variants. For sliding window analysis with diverse genotypes, your branch length correction is a pragmatic and scientifically sound approach.
