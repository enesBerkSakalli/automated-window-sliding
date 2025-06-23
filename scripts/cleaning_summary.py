#!/usr/bin/env python3
# Generate summary of alignment cleaning and preparation for sliding window analysis

from Bio import SeqIO
import os


def analyze_alignment(fasta_file):
    """Analyze basic stats of a FASTA alignment."""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if not sequences:
        return None

    return {
        "num_sequences": len(sequences),
        "alignment_length": len(sequences[0].seq),
        "total_sites": len(sequences) * len(sequences[0].seq),
    }


def main():
    data_dir = "/Users/berksakalli/Projects/automated-window-sliding/data"

    print("=== ALIGNMENT CLEANING SUMMARY ===\n")

    # Original enhanced alignment
    original_file = os.path.join(data_dir, "enhanced_alignment_combined.fasta")
    original_stats = analyze_alignment(original_file)

    # Cleaned alignment
    cleaned_file = os.path.join(data_dir, "cleaned_alignment_combined.fasta")
    cleaned_stats = analyze_alignment(cleaned_file)

    print("BEFORE CLEANING (Enhanced Combined Alignment):")
    print(f"  - Sequences: {original_stats['num_sequences']}")
    print(f"  - Alignment length: {original_stats['alignment_length']} positions")
    print(f"  - Total sites: {original_stats['total_sites']:,}")

    print("\nAFTER CLEANING:")
    print(f"  - Sequences: {cleaned_stats['num_sequences']}")
    print(f"  - Alignment length: {cleaned_stats['alignment_length']} positions")
    print(f"  - Total sites: {cleaned_stats['total_sites']:,}")

    # Calculate reductions
    seq_reduction = original_stats["num_sequences"] - cleaned_stats["num_sequences"]
    length_reduction = (
        original_stats["alignment_length"] - cleaned_stats["alignment_length"]
    )
    sites_reduction = original_stats["total_sites"] - cleaned_stats["total_sites"]

    print(f"\nREDUCTION:")
    print(f"  - Removed sequences: {seq_reduction} (VP1-only problematic sequences)")
    print(f"  - Trimmed positions: {length_reduction} (gappy columns >95% gaps)")
    print(
        f"  - Total sites reduced: {sites_reduction:,} ({sites_reduction / original_stats['total_sites'] * 100:.1f}% reduction)"
    )

    print(f"\n=== DATASET READY FOR SLIDING WINDOW ANALYSIS ===")
    print(f"Clean alignment: {cleaned_file}")
    print(
        f"Final dimensions: {cleaned_stats['num_sequences']} sequences × {cleaned_stats['alignment_length']} positions"
    )

    # Window size recommendations based on final alignment length
    alignment_len = cleaned_stats["alignment_length"]
    print(f"\nRECOMMENDED WINDOW SIZES (for {alignment_len}bp alignment):")
    print(
        f"  - Small windows: 50-75bp (high resolution, {alignment_len // 75}-{alignment_len // 50} windows)"
    )
    print(
        f"  - Medium windows: 100-150bp (balanced, {alignment_len // 150}-{alignment_len // 100} windows)"
    )
    print(
        f"  - Large windows: 200-300bp (broad patterns, {alignment_len // 300}-{alignment_len // 200} windows)"
    )

    print(f"\nNEXT STEPS:")
    print(f"1. Update parameter files to use: {cleaned_file}")
    print(f"2. Run sliding window analysis with multiple window sizes")
    print(f"3. Compare results across different window sizes")

    return cleaned_stats


if __name__ == "__main__":
    main()
