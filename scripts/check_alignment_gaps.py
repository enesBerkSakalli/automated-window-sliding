#!/usr/bin/env python3
"""
Script to analyze gaps in the enhanced alignment file.
Checks for problematic gap patterns that could affect phylogenetic analysis.
"""

from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import argparse


def analyze_alignment_gaps(fasta_file):
    """Analyze gap patterns in a multiple sequence alignment."""

    # Read sequences
    sequences = []
    seq_ids = []

    print(f"Reading alignment from {fasta_file}")
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        seq_ids.append(record.id)

    if not sequences:
        print("No sequences found in the file!")
        return

    num_seqs = len(sequences)
    seq_length = len(sequences[0])

    print(f"\nAlignment Statistics:")
    print(f"Number of sequences: {num_seqs}")
    print(f"Alignment length: {seq_length}")

    # Check if all sequences have the same length
    lengths = [len(seq) for seq in sequences]
    if len(set(lengths)) > 1:
        print(f"WARNING: Sequences have different lengths: {set(lengths)}")
        return

    # Convert to numpy array for easier analysis
    alignment_array = np.array([list(seq) for seq in sequences])

    # 1. Calculate gap percentage per sequence
    print(f"\n1. Gap Analysis per Sequence:")
    print("-" * 50)

    gap_percentages = []
    for i, seq_id in enumerate(seq_ids):
        gaps = np.sum(alignment_array[i] == "-")
        gap_percent = (gaps / seq_length) * 100
        gap_percentages.append(gap_percent)
        print(f"{seq_id[:30]:<30} {gap_percent:6.2f}% gaps ({gaps}/{seq_length})")

    # 2. Calculate gap percentage per position
    print(f"\n2. Gap Analysis per Position:")
    print("-" * 50)

    position_gaps = []
    for pos in range(seq_length):
        gaps_at_pos = np.sum(alignment_array[:, pos] == "-")
        gap_percent = (gaps_at_pos / num_seqs) * 100
        position_gaps.append(gap_percent)

    # Find positions with high gap content
    high_gap_positions = [
        (i, gap_percent)
        for i, gap_percent in enumerate(position_gaps)
        if gap_percent > 50
    ]

    print(f"Positions with >50% gaps: {len(high_gap_positions)}")
    if high_gap_positions:
        print("First 10 high-gap positions:")
        for pos, gap_percent in high_gap_positions[:10]:
            print(f"  Position {pos + 1}: {gap_percent:.1f}% gaps")

    # 3. Identify gap blocks
    print(f"\n3. Gap Block Analysis:")
    print("-" * 50)

    # Find consecutive gap regions
    gap_blocks = []
    in_gap_block = False
    block_start = 0

    for pos in range(seq_length):
        if position_gaps[pos] > 80:  # Consider >80% gaps as a gap block
            if not in_gap_block:
                block_start = pos
                in_gap_block = True
        else:
            if in_gap_block:
                gap_blocks.append((block_start, pos - 1, pos - block_start))
                in_gap_block = False

    # Handle case where alignment ends in a gap block
    if in_gap_block:
        gap_blocks.append((block_start, seq_length - 1, seq_length - block_start))

    print(f"Gap blocks (>80% gaps): {len(gap_blocks)}")
    for i, (start, end, length) in enumerate(gap_blocks):
        print(f"  Block {i + 1}: positions {start + 1}-{end + 1} (length {length})")
        if i >= 10:  # Limit output
            print(f"  ... and {len(gap_blocks) - 10} more blocks")
            break

    # 4. Overall statistics
    print(f"\n4. Overall Gap Statistics:")
    print("-" * 50)

    total_positions = num_seqs * seq_length
    total_gaps = np.sum(alignment_array == "-")
    overall_gap_percent = (total_gaps / total_positions) * 100

    print(f"Total alignment positions: {total_positions:,}")
    print(f"Total gaps: {total_gaps:,}")
    print(f"Overall gap percentage: {overall_gap_percent:.2f}%")

    # Statistics
    print(f"\nGap percentage per sequence:")
    print(f"  Min: {min(gap_percentages):.2f}%")
    print(f"  Max: {max(gap_percentages):.2f}%")
    print(f"  Mean: {np.mean(gap_percentages):.2f}%")
    print(f"  Median: {np.median(gap_percentages):.2f}%")

    print(f"\nGap percentage per position:")
    print(f"  Min: {min(position_gaps):.2f}%")
    print(f"  Max: {max(position_gaps):.2f}%")
    print(f"  Mean: {np.mean(position_gaps):.2f}%")
    print(f"  Median: {np.median(position_gaps):.2f}%")

    # 5. Quality assessment
    print(f"\n5. Alignment Quality Assessment:")
    print("-" * 50)

    # Count sequences with excessive gaps
    high_gap_seqs = sum(1 for gap_pct in gap_percentages if gap_pct > 70)
    medium_gap_seqs = sum(1 for gap_pct in gap_percentages if 30 < gap_pct <= 70)
    low_gap_seqs = sum(1 for gap_pct in gap_percentages if gap_pct <= 30)

    print(f"Sequences by gap content:")
    print(f"  Low gaps (≤30%): {low_gap_seqs} sequences")
    print(f"  Medium gaps (30-70%): {medium_gap_seqs} sequences")
    print(f"  High gaps (>70%): {high_gap_seqs} sequences")

    # Count positions with excessive gaps
    high_gap_positions_count = sum(1 for gap_pct in position_gaps if gap_pct > 80)
    medium_gap_positions_count = sum(
        1 for gap_pct in position_gaps if 50 < gap_pct <= 80
    )
    informative_positions = sum(1 for gap_pct in position_gaps if gap_pct <= 50)

    print(f"\nPositions by gap content:")
    print(
        f"  Informative (≤50% gaps): {informative_positions} positions ({informative_positions / seq_length * 100:.1f}%)"
    )
    print(
        f"  Medium gaps (50-80%): {medium_gap_positions_count} positions ({medium_gap_positions_count / seq_length * 100:.1f}%)"
    )
    print(
        f"  High gaps (>80%): {high_gap_positions_count} positions ({high_gap_positions_count / seq_length * 100:.1f}%)"
    )

    # 6. Recommendations
    print(f"\n6. Recommendations:")
    print("-" * 50)

    if overall_gap_percent > 60:
        print("⚠️  WARNING: Very high overall gap content (>60%)")
        print("   Consider trimming alignment or removing gappy sequences")
    elif overall_gap_percent > 40:
        print("⚠️  CAUTION: High overall gap content (>40%)")
        print("   Consider quality filtering before phylogenetic analysis")
    else:
        print("✅ Overall gap content is reasonable")

    if high_gap_seqs > num_seqs * 0.2:
        print(f"⚠️  WARNING: {high_gap_seqs} sequences have >70% gaps")
        print("   Consider removing these sequences")

    if informative_positions < seq_length * 0.3:
        print(f"⚠️  WARNING: Only {informative_positions} informative positions (<30%)")
        print("   Consider trimming gap-rich regions")

    if len(gap_blocks) > 10:
        print(f"⚠️  WARNING: {len(gap_blocks)} large gap blocks detected")
        print("   Consider alignment trimming or re-alignment")

    # 7. Generate gap visualization data
    print(f"\n7. Saving gap analysis data...")

    # Save detailed gap data
    with open("alignment_gap_analysis.txt", "w") as f:
        f.write("Sequence Gap Analysis\n")
        f.write("=" * 50 + "\n")
        f.write(f"Sequence_ID\tGap_Count\tGap_Percentage\n")
        for i, seq_id in enumerate(seq_ids):
            gaps = np.sum(alignment_array[i] == "-")
            gap_percent = (gaps / seq_length) * 100
            f.write(f"{seq_id}\t{gaps}\t{gap_percent:.2f}\n")

        f.write(f"\nPosition Gap Analysis\n")
        f.write("=" * 50 + "\n")
        f.write(f"Position\tGap_Count\tGap_Percentage\n")
        for pos in range(seq_length):
            gaps_at_pos = np.sum(alignment_array[:, pos] == "-")
            gap_percent = (gaps_at_pos / num_seqs) * 100
            f.write(f"{pos + 1}\t{gaps_at_pos}\t{gap_percent:.2f}\n")

    print("Gap analysis saved to 'alignment_gap_analysis.txt'")

    return {
        "num_sequences": num_seqs,
        "alignment_length": seq_length,
        "overall_gap_percent": overall_gap_percent,
        "gap_percentages": gap_percentages,
        "position_gaps": position_gaps,
        "gap_blocks": gap_blocks,
        "high_gap_seqs": high_gap_seqs,
        "informative_positions": informative_positions,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Analyze gaps in multiple sequence alignment"
    )
    parser.add_argument("fasta_file", help="Input FASTA alignment file")

    args = parser.parse_args()

    try:
        results = analyze_alignment_gaps(args.fasta_file)
        print(f"\nAnalysis complete!")

    except Exception as e:
        print(f"Error analyzing alignment: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
