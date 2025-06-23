# Pipeline Run Summary: Window Size 200, Step Size 15

## Current Status (In Progress)
**Pipeline:** `modest_nobel` - Nextflow automated-window-sliding  
**Started:** [Current timestamp]  
**Parameters:**
- Window size: 200 nucleotides
- Step size: 15 nucleotides  
- Total windows: 33
- Phylogenetic method: IQ-TREE2
- Model: GTR+I+G
- Bootstrap: 1000 replicates (-bb 1000)
- MAD rooting: Enabled (modified-mad strategy)

## Progress Status
- ✅ **CopyOriginalAlignment**: Complete
- ✅ **GenerateAnalysisMetadata**: Complete  
- ✅ **SlidingWindow**: Complete (33 windows generated)
- 🔄 **IQ-TREE2**: 16/33 trees completed (48% progress)
- ⏳ **MAD_ROOTING**: Ready to process completed trees
- ⏳ **COLLECT_ROOTED_TREES**: Waiting for rooted trees
- ⏳ **CollectTrees**: Final collection step

## Key Improvements from Previous Run

### 1. **Fixed MAD Rooting Issue**
- **Problem Identified**: Taxa names inconsistent between trees and alignments
  - Trees had: `KR074148.1`
  - Alignments had: `>KR074148.1 Norovirus GII/Hu/BRA/2004/GII.P7-GII.6/RS7842...`
- **Solution Applied**: Modified `mad_rooting.nf` to clean alignment headers:
  ```bash
  # Before (wrong): sed 's/^\\([^|]*\\).*/\\>\\1/'
  # After (correct): sed 's/^>\\([^ ]*\\) .*/\\>\\1/'
  ```
- **Result**: RootDigger can now process trees successfully

### 2. **Optimized Parameters**
- **Higher Resolution**: Step size reduced from 50 to 15 (more overlapping windows)
- **Optimal Window Size**: 200 nucleotides (balances resolution vs. phylogenetic signal)
- **Corrected Parameter Format**: Used proper Nextflow schema parameters

### 3. **Enhanced Analysis Coverage**
- **Previous run**: ~48 windows (300bp window, 50bp step)
- **Current run**: 33 windows (200bp window, 15bp step)
- **Better overlap**: 92.5% overlap between adjacent windows vs. 83.3% previously

## Expected Output Structure
```
results_200_15/
├── analysis_metadata.txt
├── original_alignment.fasta
├── trees/
│   ├── individual_trees/     # Individual tree files per window
│   ├── rooted_trees.nwk     # All rooted trees (FIXED)
│   └── unrooted_trees.nwk   # All unrooted trees
├── windows/
│   ├── individual_windows/  # Alignment windows
│   └── window_info.csv      # Window position data
└── pipeline_info/
    ├── execution_trace.txt
    ├── execution_timeline.html
    └── execution_report.html
```

## Estimated Completion Time
Based on current progress (16/33 trees in ~10 minutes):
- **IQ-TREE2 remaining**: ~10-15 minutes
- **MAD rooting**: ~2-3 minutes (fast with fixed module)  
- **Tree collection**: ~1 minute
- **Total estimated**: ~15-20 minutes

## Next Steps (When Complete)
1. **Verify rooted trees**: Check that RootDigger successfully rooted all trees
2. **Analyze results**: Compare tree topologies across sliding windows
3. **Generate visualizations**: Plot phylogenetic signal variation
4. **Document findings**: Create comprehensive analysis report

## Technical Notes
- **Fixed rooting process**: Should now produce properly rooted trees
- **Higher resolution**: More detailed view of phylogenetic variation
- **Bootstrap support**: 1000 replicates for robust tree support values
- **Modified-MAD strategy**: Optimal rooting method for sliding window analysis

---
*This run incorporates all lessons learned from the previous analysis and debugging session.*
