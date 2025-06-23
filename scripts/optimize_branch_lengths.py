#!/usr/bin/env python3
"""
Create optimized parameter file for cleaned alignment and
investigate outgroup branch length issues.
"""

import json
import os
from Bio import SeqIO


def create_cleaned_params():
    """Create parameter file for cleaned alignment."""

    # Base parameters (similar to previous successful runs)
    params = {
        "input": "/Users/berksakalli/Projects/automated-window-sliding/cleaned_alignment_with_outgroups.fasta",
        "outdir": "/Users/berksakalli/Projects/automated-window-sliding/results_cleaned_outgroups",
        "window_size": 200,
        "step_size": 15,
        "iqtree_model": "GTR+F+I+G4",
        "threads": 4,
        "bootstrap": 1000,
        "rooting_method": "mad",
    }

    # Save parameter file
    params_file = "/Users/berksakalli/Projects/automated-window-sliding/params_cleaned_outgroups.json"
    with open(params_file, "w") as f:
        json.dump(params, f, indent=2)

    print(f"Created parameter file: {params_file}")
    return params_file


def analyze_outgroup_strategy():
    """Analyze current outgroup strategy and suggest improvements."""

    alignment_file = "/Users/berksakalli/Projects/automated-window-sliding/cleaned_alignment_with_outgroups.fasta"

    if not os.path.exists(alignment_file):
        print("Cleaned alignment not found")
        return

    # Read sequences and identify outgroups
    sequences = {}
    gi_outgroups = []
    gii_sequences = []

    for record in SeqIO.parse(alignment_file, "fasta"):
        sequences[record.id] = str(record.seq)

        # Identify GI outgroups
        gi_ids = ["M87661", "L07418", "AF093797"]
        if any(gi_id in record.id for gi_id in gi_ids):
            gi_outgroups.append(record.id)
        else:
            gii_sequences.append(record.id)

    print(f"📊 OUTGROUP ANALYSIS:")
    print(f"   Total sequences: {len(sequences)}")
    print(f"   GI outgroups: {len(gi_outgroups)} ({gi_outgroups})")
    print(f"   GII ingroup: {len(gii_sequences)}")

    # Check if we should reduce outgroups to avoid long branches
    print(f"\n🎯 OUTGROUP STRATEGY OPTIONS:")
    print(f"   Current: 3 GI outgroups (M87661, L07418, AF093797)")
    print(
        f"   Option 1: Keep all 3 GI outgroups (more robust but potentially long branches)"
    )
    print(f"   Option 2: Use only 1-2 closest GI outgroups (shorter branches)")
    print(f"   Option 3: Use only M87661 (single, well-characterized outgroup)")

    return gi_outgroups, gii_sequences


def create_single_outgroup_alignment():
    """Create alignment with only one outgroup (M87661) to reduce branch lengths."""

    input_file = "/Users/berksakalli/Projects/automated-window-sliding/enhanced_alignment_with_outgroups.fasta"
    output_file = "/Users/berksakalli/Projects/automated-window-sliding/single_outgroup_alignment.fasta"

    if not os.path.exists(input_file):
        print("Input alignment not found")
        return

    kept_sequences = []
    removed_outgroups = []

    for record in SeqIO.parse(input_file, "fasta"):
        # Keep M87661 as outgroup, keep all GII sequences, remove other outgroups and KR074191.1
        if (
            "KR074191.1" in record.id
            or "L07418" in record.id
            or "AF093797" in record.id
        ):
            removed_outgroups.append(record.id)
        else:
            kept_sequences.append(record)

    # Write single outgroup alignment
    SeqIO.write(kept_sequences, output_file, "fasta")

    print(f"Created single outgroup alignment: {output_file}")
    print(f"   Kept: {len(kept_sequences)} sequences")
    print(f"   Removed: {removed_outgroups}")

    return output_file


def create_single_outgroup_params():
    """Create parameter file for single outgroup alignment."""

    params = {
        "input": "/Users/berksakalli/Projects/automated-window-sliding/single_outgroup_alignment.fasta",
        "outdir": "/Users/berksakalli/Projects/automated-window-sliding/results_single_outgroup",
        "window_size": 200,
        "step_size": 15,
        "iqtree_model": "GTR+F+I+G4",
        "threads": 4,
        "bootstrap": 1000,
        "rooting_method": "mad",
    }

    params_file = "/Users/berksakalli/Projects/automated-window-sliding/params_single_outgroup.json"
    with open(params_file, "w") as f:
        json.dump(params, f, indent=2)

    print(f"Created single outgroup parameter file: {params_file}")
    return params_file


def main():
    """Main function."""
    print("=== OPTIMIZING ALIGNMENT FOR BRANCH LENGTHS ===\n")

    # Create parameter file for cleaned alignment (without KR074191.1)
    print(
        "1. Creating parameter file for cleaned alignment (3 outgroups, no KR074191.1):"
    )
    cleaned_params = create_cleaned_params()

    # Analyze outgroup strategy
    print(f"\n2. Analyzing outgroup strategy:")
    gi_outgroups, gii_sequences = analyze_outgroup_strategy()

    # Create single outgroup alignment
    print(f"\n3. Creating single outgroup alignment:")
    single_outgroup_file = create_single_outgroup_alignment()

    # Create parameter file for single outgroup
    print(f"\n4. Creating parameter file for single outgroup:")
    single_params = create_single_outgroup_params()

    print(f"\n=== RECOMMENDATIONS ===")
    print(f"📝 Two approaches to test:")
    print(f"")
    print(f"   APPROACH 1: Cleaned 3-outgroup alignment")
    print(f"   - File: cleaned_alignment_with_outgroups.fasta (40 sequences)")
    print(f"   - Params: params_cleaned_outgroups.json")
    print(f"   - Pros: More robust rooting with multiple outgroups")
    print(f"   - Cons: May still have some long branches")
    print(f"")
    print(f"   APPROACH 2: Single outgroup alignment")
    print(f"   - File: single_outgroup_alignment.fasta (39 sequences)")
    print(f"   - Params: params_single_outgroup.json")
    print(f"   - Pros: Likely shorter, more reasonable branch lengths")
    print(f"   - Cons: Less robust rooting (single outgroup)")

    print(f"\n🚀 NEXT STEPS:")
    print(f"1. Test APPROACH 2 first (single outgroup):")
    print(
        f"   nextflow run main.nf -params-file params_single_outgroup.json -profile local"
    )
    print(f"")
    print(f"2. If branch lengths are reasonable, use single outgroup approach")
    print(f"3. If you need multiple outgroups, try APPROACH 1")
    print(f"4. Monitor branch lengths in resulting trees")


if __name__ == "__main__":
    main()
