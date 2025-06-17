#!/bin/bash

# Example: Running pipeline with exhaustive mode for higher accuracy
# Trade-off: More accurate but significantly slower

echo "Running pipeline with RootDigger exhaustive mode..."

nextflow run main.nf \
    --input "data/aligned_norovirus_sequences.fasta" \
    --outdir "test-exhaustive-accuracy" \
    --window_size 400 \
    --step_size 50 \
    --rootdigger_strategy "modified-mad" \
    --rootdigger_exhaustive true \
    --mad_rooting true \
    -work-dir "logs/work-exhaustive" \
    -resume

echo "Pipeline completed with exhaustive rooting for maximum accuracy."
echo "Check test-exhaustive-accuracy/best_rooted_trees.newick for results."
