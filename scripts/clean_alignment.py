#!/usr/bin/env python3
# Script to create cleaned alignment without problematic sequences

from Bio import SeqIO
import sys


def clean_alignment(input_file, output_file, sequences_to_remove):
    """Remove problematic sequences from alignment."""

    # Read alignment and filter out problematic sequences
    cleaned_sequences = []
    for record in SeqIO.parse(input_file, "fasta"):
        if record.id not in sequences_to_remove:
            cleaned_sequences.append(record)

    # Find columns that are mostly gaps (>95% gaps)
    if cleaned_sequences:
        alignment_length = len(cleaned_sequences[0].seq)
        columns_to_keep = []

        for i in range(alignment_length):
            gap_count = sum(1 for seq in cleaned_sequences if seq.seq[i] == "-")
            gap_percentage = gap_count / len(cleaned_sequences)

            if gap_percentage < 0.95:  # Keep columns with <95% gaps
                columns_to_keep.append(i)

        # Trim alignment to remove gappy columns
        trimmed_sequences = []
        for record in cleaned_sequences:
            trimmed_seq = "".join(record.seq[i] for i in columns_to_keep)
            record.seq = record.seq.__class__(trimmed_seq)
            trimmed_sequences.append(record)

        # Write cleaned alignment
        SeqIO.write(trimmed_sequences, output_file, "fasta")

        print(f"Cleaned alignment written to: {output_file}")
        print(f"Removed {len(sequences_to_remove)} problematic sequences")
        print(f"Trimmed from {alignment_length} to {len(columns_to_keep)} positions")
        print(
            f"Final alignment: {len(trimmed_sequences)} sequences x {len(columns_to_keep)} positions"
        )


# File paths
input_file = "/Users/berksakalli/Projects/automated-window-sliding/data/enhanced_alignment_combined.fasta"
output_file = "/Users/berksakalli/Projects/automated-window-sliding/data/cleaned_alignment_combined.fasta"

# VP1-only sequences that are problematic (have excessive leading gaps)
sequences_to_remove = {"MF681695.1", "MF681696.1"}

clean_alignment(input_file, output_file, sequences_to_remove)
