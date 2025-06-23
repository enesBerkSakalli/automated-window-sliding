# Backbone Logic Usage Guide

## Overview

The backbone logic in this pipeline provides an alternative rooting method that uses a **backbone-constrained midpoint rooting** approach. This method is faster than the default MAD rooting and provides good results for comparative phylogenetic analysis.

## How Backbone Logic Works

### 1. **Backbone Tree Generation**
- Creates a reference backbone tree from the full input alignment
- Uses IQ-TREE with automatic model selection
- Applies bootstrap support (1000 replicates)
- Applies midpoint rooting to the backbone tree

### 2. **Backbone-Constrained Window Trees**
- For each sliding window:
  - Reconstructs phylogenetic tree using IQ-TREE
  - Applies backbone constraint during tree reconstruction
  - Uses the backbone topology as guidance

### 3. **Midpoint Rooting**
- Applies fast midpoint rooting to each window tree
- Uses branch lengths to determine optimal root placement
- Much faster than exhaustive rooting methods

## Enabling Backbone Logic

### Method 1: Configuration Parameter
Edit `nextflow.config`:
```nextflow
params {
    backbone_midpoint_rooting = true    // Enable backbone logic
    mad_rooting = false                 // Disable MAD rooting (optional)
}
```

### Method 2: Command Line
```bash
nextflow run workflows/automated-window-sliding.nf \\
    --input ./data/master_alignment.fasta \\
    --backbone_midpoint_rooting true \\
    --mad_rooting false
```

### Method 3: Parameter File
Create a JSON parameter file:
```json
{
    "input": "data/master_alignment.fasta",
    "backbone_midpoint_rooting": true,
    "mad_rooting": false,
    "window_size": 300,
    "step_size": 50,
    "phylo_method": "iqtree2",
    "model": "GTR+F+I+G4"
}
```

Run with:
```bash
nextflow run workflows/automated-window-sliding.nf -params-file params_backbone.json
```

## Comparison: MAD vs Backbone Rooting

| Feature        | MAD Rooting         | Backbone Rooting    |
| -------------- | ------------------- | ------------------- |
| **Speed**      | Slower (exhaustive) | Faster (midpoint)   |
| **Accuracy**   | Maximum likelihood  | Branch length based |
| **Method**     | RootDigger          | IQ-TREE + Midpoint  |
| **Best For**   | Final analysis      | Large-scale studies |
| **Constraint** | None                | Backbone topology   |
| **Bootstrap**  | Via RootDigger      | Via backbone tree   |

## Output Files (Backbone Mode)

When backbone rooting is enabled, you'll get these additional outputs:

```
results/
├── backbone_analysis/
│   ├── backbone_tree.treefile           # Reference backbone tree
│   └── backbone_tree.log               # Generation log
├── backbone_midpoint/                   # Rooted window trees
│   ├── 1_rooted.treefile
│   ├── 1_constraint.treefile
│   ├── 101_rooted.treefile
│   └── ...
├── ALL_BACKBONE_MIDPOINT_TREES.tre     # Combined rooted trees
└── ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre
```

## Example Commands

### Basic Backbone Analysis
```bash
nextflow run workflows/automated-window-sliding.nf \\
    --input ./data/master_alignment.fasta \\
    --backbone_midpoint_rooting true
```

### Backbone with Custom Parameters
```bash
nextflow run workflows/automated-window-sliding.nf \\
    --input ./data/master_alignment.fasta \\
    --backbone_midpoint_rooting true \\
    --window_size 300 \\
    --step_size 50 \\
    --model "GTR+F+I+G4"
```

### Both MAD and Backbone (Comparison)
```bash
nextflow run workflows/automated-window-sliding.nf \\
    --input ./data/master_alignment.fasta \\
    --backbone_midpoint_rooting true \\
    --mad_rooting true
```

## When to Use Backbone Logic

### ✅ **Use Backbone Rooting When:**
- You have a large number of windows (>50)
- You need faster execution times
- You're doing preliminary analysis
- You want consistent rooting across windows
- Your sequences have reliable branch lengths

### ⚠️  **Use MAD Rooting When:**
- You need maximum accuracy
- You're doing final/publication analysis
- You have computational resources for longer runs
- Root placement is critical for your analysis

## Technical Details

### Backbone Tree Generation Process:
1. Clean alignment headers for compatibility
2. Run IQ-TREE with Model Finder Plus (MFP)
3. Generate bootstrap support (1000 replicates)
4. Apply midpoint rooting using DendroPy
5. Output rooted backbone tree

### Window Tree Processing:
1. For each window alignment:
2. Run IQ-TREE with backbone constraint
3. Apply midpoint rooting to resulting tree
4. Generate both rooted and constraint trees
5. Collect logs for quality assessment

### Dependencies:
- IQ-TREE2 (tree reconstruction)
- DendroPy (Python library for midpoint rooting)
- BioPython (sequence handling)

## Troubleshooting

### Common Issues:
1. **Missing DendroPy**: Install with `pip install dendropy`
2. **IQ-TREE errors**: Check input alignment format
3. **Memory issues**: Reduce window size or increase resources

### Performance Tips:
- Use fewer CPU cores for smaller datasets
- Increase step size for faster execution
- Use simpler evolutionary models for speed
