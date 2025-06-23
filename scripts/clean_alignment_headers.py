#!/usr/bin/env python3
"""
Clean alignment headers to contain only accession numbers.
This fixes the taxa name mismatch between alignments and trees for RootDigger.
"""

import os
import re
from Bio import SeqIO


def clean_fasta_headers(input_file, output_file):
    """
    Clean FASTA headers to contain only accession numbers.

    Args:
        input_file: Path to input FASTA file
        output_file: Path to output FASTA file with cleaned headers
    """
    cleaned_sequences = []

    with open(input_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            # Extract accession number from header (e.g., "KR074148.1" from full header)
            accession_match = re.match(r"^([A-Z]+\d+\.\d+)", record.id)
            if accession_match:
                clean_id = accession_match.group(1)
                record.id = clean_id
                record.description = ""  # Remove description to keep only accession
                cleaned_sequences.append(record)
            else:
                print(f"Warning: Could not extract accession from {record.id}")

    # Write cleaned sequences
    with open(output_file, "w") as f:
        SeqIO.write(cleaned_sequences, f, "fasta")

    print(f"Cleaned {len(cleaned_sequences)} sequences")
    print(f"Output written to: {output_file}")


def clean_all_window_alignments(results_dir):
    """
    Clean headers for all window alignment files.

    Args:
        results_dir: Path to results directory containing sliding_windows subdirectory
    """
    alignments_dir = os.path.join(results_dir, "sliding_windows", "alignments")

    if not os.path.exists(alignments_dir):
        print(f"Error: {alignments_dir} not found")
        return

    # Find all alignment files
    alignment_files = []
    for file in os.listdir(alignments_dir):
        if file.endswith(".fasta"):
            alignment_files.append(file)

    alignment_files.sort()
    print(f"Found {len(alignment_files)} alignment files to clean")

    # Create cleaned alignments directory
    cleaned_dir = os.path.join(results_dir, "sliding_windows", "alignments_cleaned")
    os.makedirs(cleaned_dir, exist_ok=True)

    # Clean each alignment file
    for alignment_file in alignment_files:
        input_path = os.path.join(alignments_dir, alignment_file)
        output_path = os.path.join(cleaned_dir, alignment_file)

        print(f"\nCleaning {alignment_file}...")
        clean_fasta_headers(input_path, output_path)


if __name__ == "__main__":
    # Clean the main refined alignment
    print("Cleaning main refined alignment...")
    clean_fasta_headers("refined_alignment.fasta", "refined_alignment_cleaned.fasta")

    # Clean all window alignments
    print("\nCleaning window alignments...")
    clean_all_window_alignments("results_optimized_200_20")

    print("\nAll alignment headers cleaned successfully!")
