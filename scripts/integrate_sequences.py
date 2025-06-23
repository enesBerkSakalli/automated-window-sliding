#!/usr/bin/env python3
"""
Script to integrate the newly fetched paper sequences with the existing norovirus sequences
and create a combined FASTA file with proper genogroup annotation.
"""

import os
import re
from Bio import SeqIO


def extract_genogroup_from_header(header):
    """Extract genogroup information from FASTA header."""
    # Look for patterns like GII.P16-GII.4 or GII.P7-GII.6
    genogroup_match = re.search(
        r"(GI{1,2}\.P[^-\s/]+(?:-GI{1,2}\.\d+|17|12|15))", header
    )
    if genogroup_match:
        return genogroup_match.group(1)

    # Look for patterns like GII.4 or GI.3
    simple_match = re.search(r"(GI{1,2}\.\d+)", header)
    if simple_match:
        return simple_match.group(1)

    return None


def standardize_header(record):
    """Standardize FASTA header format."""
    original_header = record.description
    accession = record.id

    # Extract genogroup
    genogroup = extract_genogroup_from_header(original_header)

    if genogroup:
        # Standardize the header format to match existing ones
        new_description = f"{accession} Norovirus {genogroup} - {original_header.split(' ', 1)[1] if ' ' in original_header else original_header}"
        record.description = new_description
        record.id = accession

    return record


def combine_sequences(existing_file, new_file, output_file):
    """Combine existing and new sequences into a single file."""

    print(f"Reading existing sequences from {existing_file}...")
    existing_sequences = list(SeqIO.parse(existing_file, "fasta"))
    print(f"Found {len(existing_sequences)} existing sequences")

    print(f"Reading new sequences from {new_file}...")
    new_sequences = list(SeqIO.parse(new_file, "fasta"))
    print(f"Found {len(new_sequences)} new sequences")

    # Standardize new sequence headers
    print("Standardizing new sequence headers...")
    standardized_new_sequences = []
    for record in new_sequences:
        standardized_record = standardize_header(record)
        standardized_new_sequences.append(standardized_record)

        # Print genogroup info for first few sequences
        if len(standardized_new_sequences) <= 5:
            genogroup = extract_genogroup_from_header(record.description)
            print(f"  {record.id}: {genogroup if genogroup else 'No genogroup found'}")

    # Combine all sequences
    all_sequences = existing_sequences + standardized_new_sequences

    print(f"Writing {len(all_sequences)} total sequences to {output_file}...")
    with open(output_file, "w") as f:
        SeqIO.write(all_sequences, f, "fasta")

    return len(existing_sequences), len(standardized_new_sequences), len(all_sequences)


def create_genogroup_summary(combined_file):
    """Create a summary of genogroups in the combined file."""
    genogroup_counts = {}
    accession_to_genogroup = {}

    for record in SeqIO.parse(combined_file, "fasta"):
        genogroup = extract_genogroup_from_header(record.description)
        accession = record.id

        if genogroup:
            genogroup_counts[genogroup] = genogroup_counts.get(genogroup, 0) + 1
            accession_to_genogroup[accession] = genogroup
        else:
            print(f"Warning: No genogroup found for {accession}")

    return genogroup_counts, accession_to_genogroup


def main():
    """Main function to integrate sequences."""

    # Configuration
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    existing_file = os.path.join(base_dir, "data", "norovirus_sequences.fasta")
    new_file = os.path.join(base_dir, "data", "paper_sequences.fasta")
    combined_file = os.path.join(base_dir, "data", "combined_norovirus_sequences.fasta")
    summary_file = os.path.join(base_dir, "data", "genogroup_summary.txt")

    # Check if files exist
    if not os.path.exists(existing_file):
        print(f"Error: Existing file not found: {existing_file}")
        return

    if not os.path.exists(new_file):
        print(f"Error: New sequences file not found: {new_file}")
        return

    # Combine sequences
    existing_count, new_count, total_count = combine_sequences(
        existing_file, new_file, combined_file
    )

    # Create genogroup summary
    print("Analyzing genogroups in combined file...")
    genogroup_counts, accession_to_genogroup = create_genogroup_summary(combined_file)

    # Write summary to file
    with open(summary_file, "w") as f:
        f.write("Genogroup Summary for Combined Norovirus Sequences\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total sequences: {total_count}\n")
        f.write(f"  - Existing sequences: {existing_count}\n")
        f.write(f"  - New sequences from paper: {new_count}\n\n")

        f.write("Genogroup distribution:\n")
        for genogroup, count in sorted(genogroup_counts.items()):
            f.write(f"  {genogroup}: {count} sequences\n")

        f.write(f"\nTotal unique genogroups: {len(genogroup_counts)}\n")

    # Print summary
    print("\nIntegration Summary:")
    print(f"- Existing sequences: {existing_count}")
    print(f"- New sequences from paper: {new_count}")
    print(f"- Total combined sequences: {total_count}")
    print(f"- Combined file saved to: {combined_file}")
    print(f"- Summary saved to: {summary_file}")

    print("\nGenogroup distribution:")
    for genogroup, count in sorted(genogroup_counts.items()):
        print(f"  {genogroup}: {count} sequences")

    print(f"\nTotal unique genogroups: {len(genogroup_counts)}")


if __name__ == "__main__":
    main()
