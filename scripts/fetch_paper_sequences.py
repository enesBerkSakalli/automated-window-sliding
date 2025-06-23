#!/usr/bin/env python3
"""
Script to identify missing accession numbers from the PLoS ONE paper (10.1371/journal.pone.0189504)
and fetch them from GenBank to integrate with the existing norovirus sequences.
"""

import os
import re
from Bio import Entrez, SeqIO
import time


def get_paper_accession_numbers():
    """Get all accession numbers mentioned in the PLoS ONE paper."""
    accession_ranges = [
        ("KY451971", "KY451987"),  # 17 sequences
        ("MF158177", "MF158199"),  # 23 sequences
        ("MF681695", "MF681696"),  # 2 sequences
        ("KY551568", "KY551569"),  # 2 sequences
    ]

    all_accessions = []

    for start, end in accession_ranges:
        # Extract the numeric part and prefix
        start_num = int(start[2:])  # Remove the 2-letter prefix
        end_num = int(end[2:])
        prefix = start[:2]

        for num in range(start_num, end_num + 1):
            accession = f"{prefix}{num:06d}.1"  # Add .1 version
            all_accessions.append(accession)

    return sorted(all_accessions)


def get_current_accessions(fasta_file):
    """Extract accession numbers from the current FASTA file."""
    current_accessions = set()

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Extract accession number (first part of header)
                match = re.match(r">(\w+\.\d+)", line)
                if match:
                    current_accessions.add(match.group(1))

    return current_accessions


def fetch_sequences_from_genbank(accession_list, email):
    """Fetch sequences from GenBank for the given accession numbers."""
    Entrez.email = email
    sequences = []

    print(f"Fetching {len(accession_list)} sequences from GenBank...")

    for i, accession in enumerate(accession_list):
        try:
            print(f"Fetching {accession} ({i + 1}/{len(accession_list)})...")

            # Fetch the sequence
            handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="fasta", retmode="text"
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()

            sequences.append(record)

            # Be nice to NCBI servers
            time.sleep(0.5)

        except Exception as e:
            print(f"Error fetching {accession}: {e}")
            continue

    return sequences


def main():
    """Main function to identify and fetch missing sequences."""

    # Configuration
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    current_fasta_file = os.path.join(base_dir, "data", "norovirus_sequences.fasta")
    output_file = os.path.join(base_dir, "data", "paper_sequences.fasta")

    print("Getting accession numbers from the PLoS ONE paper...")
    paper_accessions = get_paper_accession_numbers()
    print(f"Found {len(paper_accessions)} accession numbers in the paper")

    print("Getting current accession numbers from FASTA file...")
    current_accessions = get_current_accessions(current_fasta_file)
    print(f"Found {len(current_accessions)} sequences in current file")

    # Find missing accessions
    missing_accessions = [
        acc for acc in paper_accessions if acc not in current_accessions
    ]

    if not missing_accessions:
        print("All sequences from the paper are already present in your FASTA file!")
        return

    print(f"\nMissing accession numbers ({len(missing_accessions)}):")
    for acc in missing_accessions:
        print(f"  {acc}")

    # Ask user for email before proceeding
    user_email = input("\nPlease enter your email address for NCBI access: ").strip()
    if not user_email or "@" not in user_email:
        print("Invalid email address. Please run the script again with a valid email.")
        return

    # Fetch missing sequences
    print(f"\nFetching {len(missing_accessions)} missing sequences from GenBank...")
    missing_sequences = fetch_sequences_from_genbank(missing_accessions, user_email)

    if missing_sequences:
        # Write sequences to file
        with open(output_file, "w") as f:
            SeqIO.write(missing_sequences, f, "fasta")

        print(f"\nSuccessfully fetched {len(missing_sequences)} sequences!")
        print(f"Sequences saved to: {output_file}")

        # Show summary
        print("\nSummary:")
        print(f"- Paper accessions: {len(paper_accessions)}")
        print(f"- Current file: {len(current_accessions)}")
        print(f"- Missing: {len(missing_accessions)}")
        print(f"- Successfully fetched: {len(missing_sequences)}")

        if len(missing_sequences) < len(missing_accessions):
            failed = len(missing_accessions) - len(missing_sequences)
            print(f"- Failed to fetch: {failed}")

    else:
        print("No sequences were successfully fetched.")


if __name__ == "__main__":
    main()
