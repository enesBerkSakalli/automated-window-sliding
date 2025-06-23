# Results Directory Structure

## Overview

The Automated Window Sliding pipeline automatically organizes all results in a structured directory hierarchy under the main `results/` folder. This ensures consistent organization, reproducibility, and easy access to analysis outputs.

## Directory Structure

```
results/
├── [INPUT_NAME]_w[WINDOW_SIZE]_s[STEP_SIZE]_[MODEL]_[TIMESTAMP]/
│   ├── analysis_metadata.txt                 # Run parameters and metadata
│   ├── original_alignment_[INPUT_NAME].fasta # Copy of input alignment
│   ├── windows.log                           # Sliding window information
│   ├── tree_reconstruction_logs/             # Detailed logs from tree reconstruction
│   │   ├── iqtree_files/                     # IQ-TREE output files
│   │   └── logs/                            # Process logs
│   └── rooted_trees/                        # MAD rooted trees (if enabled)
├── best_trees.newick                        # Combined best trees (Newick format)
├── best_trees.nexus                         # Combined best trees (Nexus format)
├── best_rooted_trees.newick                 # Combined rooted trees (if MAD rooting enabled)
├── best_rooted_trees.nexus                  # Combined rooted trees (if MAD rooting enabled)
└── pipeline_info/                           # Pipeline execution reports
    ├── execution_timeline_[TIMESTAMP].html  # Timeline visualization
    ├── execution_report_[TIMESTAMP].html    # Execution report
    ├── execution_trace_[TIMESTAMP].txt      # Detailed trace
    └── pipeline_dag_[TIMESTAMP].html        # DAG visualization
```

## Key Features

### 1. Automatic Organization
- **Main results folder**: All outputs go to `results/` by default
- **Timestamped runs**: Each analysis creates a unique timestamped subdirectory
- **Structured hierarchy**: Organized by analysis type and output format

### 2. Reproducibility
- **Input preservation**: Original alignment file copied to results
- **Metadata tracking**: Complete analysis parameters recorded
- **Version control**: Timestamped directories prevent overwrites

### 3. Easy Access
- **Summary files**: Combined results at top level for quick access
- **Detailed files**: Individual window results in timestamped subdirectories
- **Multiple formats**: Newick and Nexus formats provided

## Configuration

The results directory structure is controlled by these parameters:

```bash
# Default results directory (automatically set)
--outdir results

# Custom results directory
--outdir custom_results_folder
```

## Example Usage

```bash
# Standard run - results go to results/
nextflow run main.nf \\
  --input alignment.fasta \\
  --window_size 200 \\
  --step_size 50 \\
  --model "GTR+F+I+G4"

# Custom output directory - results go to my_analysis/
nextflow run main.nf \\
  --input alignment.fasta \\
  --outdir my_analysis \\
  --window_size 200 \\
  --step_size 50
```

## Output Files Description

### Summary Files (in main results/ directory)
- **`best_trees.newick/nexus`**: All reconstructed trees combined
- **`best_rooted_trees.newick/nexus`**: All MAD rooted trees (if enabled)

### Run-Specific Directory
- **`analysis_metadata.txt`**: Complete run parameters and system information
- **`original_alignment_*.fasta`**: Input alignment for reproducibility
- **`windows.log`**: Sliding window positions and statistics
- **`tree_reconstruction_logs/`**: Detailed IQ-TREE/RAxML output files

### Pipeline Information
- **`pipeline_info/`**: Nextflow execution reports and visualizations

## Tips

1. **Multiple runs**: Each run creates a new timestamped directory
2. **Comparison**: Use timestamps to compare different parameter sets
3. **Cleanup**: Remove old timestamped directories to save space
4. **Backup**: The structured format makes backup and archiving easy

---

*This organization ensures your sliding window phylogenetic analyses are reproducible, well-documented, and easy to navigate.*
