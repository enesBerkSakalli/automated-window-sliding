#!/usr/bin/env python3
"""
Script to fix alignment files by truncating sequence names to match tree format.
This creates 'fixed' versions of alignment files with shortened sequence names.
"""

import os
import sys
from Bio import SeqIO
from pathlib import Path

def fix_alignment_names(input_file, output_file):
    """Fix sequence names in alignment file to match tree format."""
    
    fixed_sequences = []
    
    with open(input_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Extract only the accession number (before first space)
            new_id = record.id.split()[0]
            record.id = new_id
            record.description = ""  # Remove description
            fixed_sequences.append(record)
    
    with open(output_file, 'w') as f:
        SeqIO.write(fixed_sequences, f, 'fasta')
    
    print(f"Fixed {len(fixed_sequences)} sequences: {input_file} -> {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python fix_alignment_names.py <input_alignment> <output_alignment>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)
    
    fix_alignment_names(input_file, output_file)

if __name__ == "__main__":
    main()
