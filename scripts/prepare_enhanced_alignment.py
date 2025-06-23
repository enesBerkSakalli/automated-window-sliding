#!/usr/bin/env python3
"""
Prepare enhanced alignment from combined norovirus sequences for sliding window analysis
"""

import os
import subprocess
from Bio import SeqIO


def create_enhanced_alignment():
    """Create alignment from combined sequences using MAFFT."""

    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    combined_fasta = os.path.join(
        base_dir, "data", "combined_norovirus_sequences.fasta"
    )
    output_alignment = os.path.join(
        base_dir, "data", "enhanced_alignment_combined.fasta"
    )

    print("Creating enhanced alignment from combined sequences...")
    print(f"Input: {combined_fasta}")
    print(f"Output: {output_alignment}")

    # Check if input exists
    if not os.path.exists(combined_fasta):
        print(f"Error: Combined sequences file not found: {combined_fasta}")
        return False

    # Count sequences
    seq_count = len(list(SeqIO.parse(combined_fasta, "fasta")))
    print(f"Sequences to align: {seq_count}")

    # Run MAFFT alignment
    try:
        print("Running MAFFT alignment (this may take a few minutes)...")
        cmd = ["mafft", "--auto", "--thread", "4", combined_fasta]

        with open(output_alignment, "w") as outfile:
            result = subprocess.run(
                cmd, stdout=outfile, stderr=subprocess.PIPE, text=True
            )

        if result.returncode == 0:
            print("✅ Alignment completed successfully!")
            print(f"Output saved to: {output_alignment}")

            # Check output
            aligned_seqs = list(SeqIO.parse(output_alignment, "fasta"))
            print(f"Aligned sequences: {len(aligned_seqs)}")
            if aligned_seqs:
                print(f"Alignment length: {len(aligned_seqs[0].seq)} bp")

            return True
        else:
            print("❌ MAFFT failed with error:")
            print(result.stderr)
            return False

    except FileNotFoundError:
        print("❌ MAFFT not found. Please install MAFFT first:")
        print("  conda install -c bioconda mafft")
        print("  or")
        print("  brew install mafft")
        return False
    except Exception as e:
        print(f"❌ Error running MAFFT: {e}")
        return False


def main():
    """Main function to create enhanced alignment."""
    print("🔬 PREPARING ENHANCED ALIGNMENT FOR SLIDING WINDOW ANALYSIS")
    print("=" * 60)

    success = create_enhanced_alignment()

    if success:
        print("\n🎉 Enhanced alignment is ready for sliding window analysis!")
        print("\nNext steps:")
        print("1. Parameter files will be created for three window sizes")
        print("2. Sliding window analysis will be run with enhanced dataset")
        print("3. Results will include 88 sequences with 17 genogroups")
    else:
        print(
            "\n❌ Alignment preparation failed. Please check the error messages above."
        )


if __name__ == "__main__":
    main()
