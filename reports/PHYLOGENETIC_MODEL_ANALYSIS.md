# Phylogenetic Model Analysis Report

## Pipeline Run Details
- **Window Size**: 200 bp
- **Step Size**: 15 bp  
- **Number of Windows**: 33
- **Run Date**: 2025-06-18
- **Results Directory**: results_200_15

## Model Configuration vs. Actual Model Used

### Parameters Specified
In `params_200_15.json`, we specified:
```json
{
    "model": "GTR+I+G",
    "phylo_method": "iqtree2",
    "phylo_parameters": "-bb 1000"
}
```

### Actual Model Used
According to the IQ-TREE log files, **all 33 sliding windows** used the model:
**GTR+F+I+G4**

### Model Components Explanation

**GTR+F+I+G4** consists of:
- **GTR**: General Time Reversible substitution model (6 rate parameters)
- **+F**: Empirical base frequencies (estimated from the data)
- **+I**: Proportion of invariable sites
- **+G4**: Gamma rate heterogeneity with 4 discrete categories

### Differences from Specified Model

| Component          | Specified     | Actual | Explanation                         |
| ------------------ | ------------- | ------ | ----------------------------------- |
| Base Model         | GTR           | GTR    | ✅ Same                              |
| Base Frequencies   | Not specified | +F     | IQ-TREE added empirical frequencies |
| Invariable Sites   | +I            | +I     | ✅ Same                              |
| Rate Heterogeneity | +G            | +G4    | IQ-TREE used 4 categories (default) |

### Why the Model Changed

1. **Automatic Model Enhancement**: IQ-TREE automatically enhances the specified model by:
   - Adding empirical base frequencies (+F) when not specified
   - Using 4 discrete gamma categories (+G4) as the default for +G

2. **Model Selection vs. Fixed Model**: Even when a model is specified, IQ-TREE may:
   - Optimize the model parameters
   - Add components that improve the likelihood
   - Use ModelFinder to validate the model choice

### Model Quality Assessment

From the IQ-TREE logs (example from window 271):
- **Log-likelihood**: -1605.6291
- **AIC score**: 3357.2582
- **BIC score**: 3598.0354
- **Free parameters**: 73 (branches + model parameters)

### Rate Parameters (Example from Window 271)
```
Rate parameter R:
  A-C: 0.7267
  A-G: 2.3606
  A-T: 1.7200
  C-G: 0.8076
  C-T: 3.6986
  G-T: 1.0000

State frequencies:
  pi(A) = 0.2347
  pi(C) = 0.2576
  pi(G) = 0.293
  pi(T) = 0.2147
```

### Gamma Rate Heterogeneity
- **Shape parameter (α)**: ~0.7457 (indicating moderate rate variation)
- **Proportion of invariable sites**: ~7.95e-05 (very small, almost negligible)

## Conclusions

1. **Consistent Model Usage**: All 33 windows used the same sophisticated model (GTR+F+I+G4)

2. **Appropriate Model Choice**: The GTR+F+I+G4 model is excellent for norovirus sequences as it:
   - Accounts for all possible nucleotide substitution rates (GTR)
   - Uses empirical base frequencies from the data (+F)
   - Handles invariable sites (+I)
   - Models rate heterogeneity across sites (+G4)

3. **Enhanced from Specification**: IQ-TREE intelligently enhanced our GTR+I+G specification to GTR+F+I+G4, which is more realistic for real sequence data

4. **High-Quality Phylogenetic Inference**: The consistent use of this sophisticated model across all windows ensures reliable and comparable phylogenetic estimates

## Recommendations

For future runs, you could specify the full model explicitly:
```json
{
    "model": "GTR+F+I+G4",
    "phylo_method": "iqtree2",
    "phylo_parameters": "-bb 1000"
}
```

However, IQ-TREE's automatic enhancement from GTR+I+G to GTR+F+I+G4 demonstrates its intelligent model selection capabilities, so the current configuration is perfectly acceptable.
