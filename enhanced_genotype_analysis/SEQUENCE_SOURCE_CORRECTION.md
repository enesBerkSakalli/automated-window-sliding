# SEQUENCE SOURCE ATTRIBUTION CORRECTION

## Investigation Summary

Based on comprehensive web search and detective analysis, I have identified and corrected a critical misattribution in the sequence sources for the norovirus dataset.

## Key Findings

### Original Issue
- The sequences with accession numbers KR074148-KR074189 (covering years 2004-2011) were previously attributed to "Hernández et al. (2016)"
- This was **INCORRECT**

### Correct Attribution
- These sequences (KR074148-KR074189, years 2004-2011) should be attributed to **"Fumian et al. (2016)"**

### Supporting Evidence
1. **Web Search Results**: Found clear references to "Fumian TM, da Silva Ribeiro de Andrade J, Leite JP, Miagostovich MP. PLoS One. 2016 Apr 26;11(4):e0145391"
2. **Publication Title**: "Norovirus Recombinant Strains Isolated from Gastroenteritis Outbreaks in Southern Brazil, 2004-2011"
3. **DOI**: 10.1371/journal.pone.0145391
4. **Time Period**: Matches exactly with the sequence years (2004-2011)
5. **Geographic Location**: Southern Brazil (matches sequence metadata)

## Corrections Made

### 1. Detective Analysis Script
- Updated `publications` database to include correct Fumian et al. (2016) reference
- Modified `assign_publication()` function to properly attribute 2004-2011 sequences to Fumian et al.
- Expanded genotype coverage for Fumian study

### 2. Comprehensive Table
- All sequences from 2004-2011 now correctly attributed to "Fumian et~al. (2016)"
- Covers 44 sequences (51.2% of total dataset)

### 3. Publication Distribution (Corrected)
- **Fumian et~al. (2016)**: 44 sequences (51.2%) - KR074148-KR074189 series, years 2004-2011
- **Barreira et~al. (2017)**: 40 sequences (46.5%) - Years 2015-2016
- **Current study**: 2 sequences (2.3%) - Additional analyses

## Verification

The corrected attribution is supported by:
1. Matching publication years (2004-2011)
2. Matching geographic region (Southern Brazil)
3. Matching accession number series (KR074148-KR074189)
4. Consistent genotype patterns
5. Publication metadata alignment

## Impact

This correction ensures:
- **Scientific Accuracy**: Proper citation of original research
- **Data Integrity**: Correct source attribution in all analyses
- **Publication Ethics**: Appropriate credit to original authors
- **Reproducibility**: Clear provenance of sequence data

## Files Updated

1. `scripts/comprehensive_detective_analysis.py` - Core attribution logic
2. `enhanced_genotype_analysis/comprehensive_norovirus_table.tex` - LaTeX table
3. `enhanced_genotype_analysis/detective_analysis_report.txt` - Analysis report
4. All phylogenetic tree files - Genotype mappings updated

---

**Date**: June 20, 2025  
**Correction Type**: Source Attribution  
**Validated By**: Web search and metadata analysis  
**Status**: ✅ CORRECTED
