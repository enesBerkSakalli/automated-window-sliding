# 🧬 Genogroup Mapping Completion Report

## 🎯 Task Completed Successfully

The genogroup mapping script has been successfully applied to all tree files in the specified results directory:
`results/cleaned_alignment_combined_w250_s10_GTR_F_I_G4_20250619_103929`

---

## 📊 Processing Summary

### ✅ Files Processed
- **Main Tree Collections**: 3 files processed
  - `all_rooted_trees.newick` 
  - `all_constraint_trees.newick`
  - `rooted_trees_collection.newick`

- **Individual Window Trees**: 104 files processed
  - **52 rooted trees** (positions 1, 11, 21, ..., 511)
  - **52 constraint trees** (positions 1, 11, 21, ..., 511)

### 🧬 Genogroup Diversity
- **86 sequences** successfully mapped
- **17 unique genogroups** identified:
  ```
  GII.1, GII.3, GII.4, GII.P13-GII.17, GII.P15-GII.15, 
  GII.P16-GII.3, GII.P16-GII.4, GII.P17_GII.17, GII.P2-GII.2, 
  GII.P21-GII.13, GII.P21-GII.21, GII.P21-GII.3, GII.P4-GII.4, 
  GII.P7-GII.14, GII.P7-GII.6, GII.Pe-GII.17, GII.Pg-GII.12
  ```

---

## 📁 Output Files Created

### 1. Reference Mapping File
- **`accession_genogroup_mapping.txt`** - Complete lookup table of all 86 sequences

### 2. Main Tree Collections (3 files each format)
- **`*_with_genogroups.treefile`** - Trees with `genogroup_accession` format
- **`*_genogroups_only.treefile`** - Trees with genogroup names only

### 3. Individual Window Trees (52 windows × 2 types × 2 formats = 208 files)
- **Rooted trees**: `{window}_rooted_genogroups_only.treefile` and `{window}_rooted_with_genogroups.treefile`
- **Constraint trees**: `{window}_constraint_genogroups_only.treefile` and `{window}_constraint_with_genogroups.treefile`

---

## 🔍 Sample Mappings

| Accession  | Genogroup     | Source Study                |
| ---------- | ------------- | --------------------------- |
| KR074148.1 | GII.P7-GII.6  | Southern Brazil (2004-2011) |
| KR074149.1 | GII.P7-GII.14 | Southern Brazil (2004-2011) |
| KY451971.1 | GII.P16-GII.4 | Espírito Santo (2015-2016)  |
| MF158177.1 | GII.P16-GII.4 | Espírito Santo (2015-2016)  |

---

## 🎯 Usage Guide

### For Visualization
- **Use `*_genogroups_only.treefile`** files for clean visualization
- Example: `1_rooted_genogroups_only.treefile` shows genogroups like `GII.P7-GII.6`, `GII.P16-GII.4`

### For Detailed Analysis  
- **Use `*_with_genogroups.treefile`** files to maintain traceability
- Example: `1_rooted_with_genogroups.treefile` shows `GII.P7-GII.6_KR074148.1`

### For Reference
- **Use `accession_genogroup_mapping.txt`** for complete lookup table

---

## 📈 Analysis Applications

### 🔬 Recombination Analysis
- **Window-by-window comparison**: Use individual window trees to track genogroup clustering changes
- **Breakpoint detection**: Compare genogroup relationships across adjacent windows
- **Phylogenetic incongruence**: Identify windows where genogroup relationships differ significantly

### 🌳 Visualization
- **Clean tree displays**: Genogroup-only trees for publication figures
- **Comparative analysis**: Side-by-side window comparisons
- **Temporal dynamics**: Track genogroup evolution across the 12-year study period

### 📊 Statistical Analysis
- **Robinson-Foulds distances**: Calculate topological differences between windows
- **Support value analysis**: Assess confidence in genogroup relationships
- **Clustering patterns**: Identify consistent vs. variable genogroup associations

---

## 🎯 Next Steps Recommendations

### 1. Topological Analysis
```bash
# Example workflow for next phase
python analyze_topological_incongruence.py \
  --input-dir "results/cleaned_alignment_combined_w250_s10_GTR_F_I_G4_20250619_103929/backbone_midpoint" \
  --tree-pattern "*_rooted_genogroups_only.treefile" \
  --output "topological_analysis_results"
```

### 2. Visualization
- Load `*_genogroups_only.treefile` files into tree visualization software
- Create sliding window phylogenetic plots
- Generate heatmaps of genogroup relationships

### 3. Publication Preparation
- Use genogroup-labeled trees for manuscript figures
- Reference the mapping file for materials and methods
- Include genogroup diversity statistics in results

---

## 🏆 Success Metrics

- ✅ **100% Processing Success**: All 107 tree files successfully processed
- ✅ **Complete Coverage**: All 86 sequences mapped to genogroups
- ✅ **High Diversity**: 17 unique genogroups preserved
- ✅ **Multiple Formats**: Both detailed and simplified tree versions created
- ✅ **Full Traceability**: Complete mapping documentation provided

---

## 📋 File Locations

**Main Results Directory:**
```
results/cleaned_alignment_combined_w250_s10_GTR_F_I_G4_20250619_103929/
├── accession_genogroup_mapping.txt
├── all_rooted_trees_genogroups_only.treefile
├── all_rooted_trees_with_genogroups.treefile
├── all_constraint_trees_genogroups_only.treefile
├── all_constraint_trees_with_genogroups.treefile
├── rooted_trees_collection_genogroups_only.treefile
├── rooted_trees_collection_with_genogroups.treefile
└── backbone_midpoint/
    ├── 1_rooted_genogroups_only.treefile
    ├── 1_rooted_with_genogroups.treefile
    ├── 1_constraint_genogroups_only.treefile
    ├── 1_constraint_with_genogroups.treefile
    └── [... 51 more windows with same pattern]
```

---

**🎉 The genogroup mapping is now complete and ready for advanced phylogenetic analysis, visualization, and publication preparation!**

*Processing completed: [Current Date]*  
*Total files created: 215 (1 mapping + 214 tree files)*  
*Script used: `map_genogroups_custom.py`*
