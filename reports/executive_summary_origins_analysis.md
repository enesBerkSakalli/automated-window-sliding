# Executive Summary: Brazilian Norovirus Recombination Dataset Analysis

## Dataset Overview
- **86 high-quality sequences** × **516 bp** (ORF1/ORF2 junction)
- **17 unique genogroups** representing extensive Brazilian norovirus diversity
- **12-year temporal span** (2004-2016) covering **critical** evolutionary periods
- **Two geographic regions**: Southern Brazil and Espírito Santo State

## Key Sequence Sources

### **1**. Southern Brazil Collection (2004-2011)
- **44 sequences** from KR074148-KR074191
- **First comprehensive report** of 8 different recombinant strains in Brazil
- **Fumian et al. 2016**, PLoS One 11(4):e0145391
- Focus on gastroenteritis outbreak strains

### 2. Espírito Santo State Collection (2015-2016)
- **42 sequences** from multiple accession ranges
- **Novel GII.P16-GII.4 Sydney variant** detection
- **PLoS One 2017**, 12(12):e0189504
- Contemporary recombinant emergence

## Analysis Specifications

### Primary Analysis Parameters
- **Window Size**: 250 bp (optimal phylogenetic signal)
- **Step Size**: 10 bp (high-resolution breakpoint detection)
- **Method**: IQ-TREE2 with GTR+F+I+G4 model
- **Rooting**: Backbone midpoint for topological stability

### Quality Control Applied
- ✅ Removed 2 problematic VP1-only sequences
- ✅ Trimmed gaps: 1,857 → 516 positions
- ✅ Preserved 49.22% phylogenetically informative sites
- ✅ Maintained complete genogroup diversity

## Research Significance

### Analytical Strengths
🌟 **Recombination Hotspot Focus**: ORF1/ORF2 junction sequences  
🌟 **Temporal Dynamics**: 12-year evolutionary coverage  
🌟 **Geographic Sampling**: Multiple Brazilian circulation areas  
🌟 **Quality Assurance**: Peer-reviewed publication sources  
🌟 **Methodological Rigor**: State-of-the-art phylogenetic inference  

### Expected Outcomes
- High-resolution recombination breakpoint mapping
- Temporal dynamics of recombinant strain emergence
- Regional patterns of norovirus circulation
- Enhanced surveillance and outbreak preparedness insights

## Impact and Applications

### Scientific Contributions
- Comprehensive Brazilian norovirus recombination characterization
- Optimized sliding window analysis methodology
- Integration of sequence origins with phylogenetic analysis

### Public Health Implications
- Improved understanding of recombinant strain circulation
- Enhanced surveillance capabilities
- Vaccine development considerations

---

**This dataset represents the most comprehensive collection of Brazilian norovirus recombinant sequences available for high-resolution phylogenetic analysis, combining rigorous sequence provenance research with state-of-the-art analytical methods.**

---

*Analysis Framework: Nextflow sliding window phylogenetic pipeline*  
*Primary Output: results_cleaned_large_windows_250_10/*  
*Full Report: comprehensive_origins_analysis_report.md*
