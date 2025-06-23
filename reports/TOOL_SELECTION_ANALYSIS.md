# Tool Selection Analysis: IQ-TREE vs RAxML-ng for Norovirus Sliding Window Analysis

## 🔍 Research Summary

Based on current literature (2024-2025), here's the comparison for your specific use case:

### Dataset Characteristics
- **Sequences:** 38 norovirus sequences  
- **Length:** 484 bp (after enhancement)
- **Analysis:** Sliding window phylogenetics
- **Target:** ~20-25 windows of 75-100 bp each

### Tool Comparison

| Factor              | IQ-TREE                | RAxML-ng           |
| ------------------- | ---------------------- | ------------------ |
| **Availability**    | ✅ Installed            | ❌ Not installed    |
| **Model Selection** | ✅ Superior ModelFinder | ⚠️ Basic options    |
| **Small Datasets**  | ✅ Excellent stability  | ✅ Good performance |
| **Speed**           | ✅ Fast                 | ✅ Very fast        |
| **Viral Studies**   | ✅ Widely used          | ✅ Also used        |
| **Window Analysis** | ✅ Optimized            | ✅ Good             |

### Literature Evidence

1. **Recent Norovirus Studies (2024-2025):**
   - Use ModelFinder (IQ-TREE) for model selection
   - Both tools used for ML inference
   - IQ-TREE more common in viral phylogenetics

2. **Performance on Small Datasets:**
   - Both tools perform similarly well on datasets <100 sequences
   - IQ-TREE shows better stability
   - RAxML-ng advantage mainly on large datasets (>1000 sequences)

3. **Sliding Window Analysis:**
   - RAxML classic has built-in sliding window support
   - IQ-TREE works excellently with external window pipelines
   - Both suitable for your use case

## 🎯 Recommendation: USE IQ-TREE

### Reasons:
1. **Already installed** - No setup time needed
2. **Superior model selection** - ModelFinder is industry standard
3. **Stable on small datasets** - Perfect for your 38 sequences
4. **Proven track record** - Widely used in recent viral studies
5. **Integration ready** - Works seamlessly with the pipeline

### Pipeline Configuration:
```bash
# Use IQ-TREE for optimal results
--phylo_method iqtree2
--model_criterion bic
```

### If you want RAxML-ng later:
```bash
# Install via Homebrew
brew install brewsci/bio/raxml-ng

# Or via Conda
conda install -c bioconda raxml-ng
```

## ✅ Conclusion

For your norovirus sliding window analysis:
- **IQ-TREE is the optimal choice** given your dataset size and requirements
- The performance difference would be minimal for your 38-sequence dataset
- ModelFinder gives IQ-TREE a significant advantage
- You can always install RAxML-ng later for comparison

**Proceed with IQ-TREE for your analysis!** 🚀
