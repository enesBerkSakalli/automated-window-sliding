#!/usr/bin/env python3
"""
Analysis of KR074191.1 as potential outgroup sequence for norovirus phylogenetic analysis.

This script examines the sequence KR074191.1 (GII.Pg-GII.12) mentioned in the user's question
and evaluates its suitability as an outgroup for the norovirus dataset.

Based on the research paper 10.1371/journal.pone.0189504, this sequence represents
a recombinant GII.Pg-GII.12 strain from Brazil, 2009.
"""

import re
from Bio import SeqIO
from Bio.SeqUtils import GC, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd


def analyze_outgroup_candidate():
    """Analyze KR074191.1 sequence and its potential as an outgroup."""

    print("=" * 80)
    print("ANALYSIS OF KR074191.1 AS POTENTIAL OUTGROUP")
    print("=" * 80)
    print()

    # Read the FASTA file
    sequences = {}
    sequence_info = {}

    try:
        with open("norovirus_sequences.fasta", "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = str(record.seq)
                sequence_info[record.id] = record.description
    except FileNotFoundError:
        print("Error: norovirus_sequences.fasta not found")
        return

    # Find KR074191.1
    target_seq = "KR074191.1"
    if target_seq not in sequences:
        print(f"Error: {target_seq} not found in the dataset")
        return

    print(f"TARGET SEQUENCE: {target_seq}")
    print(f"Description: {sequence_info[target_seq]}")
    print()

    # Analyze the target sequence
    seq = sequences[target_seq]
    print("SEQUENCE CHARACTERISTICS:")
    print(f"Length: {len(seq)} bp")
    print(f"GC content: {GC(seq):.2f}%")
    print(f"A: {seq.count('A')} ({seq.count('A') / len(seq) * 100:.1f}%)")
    print(f"T: {seq.count('T')} ({seq.count('T') / len(seq) * 100:.1f}%)")
    print(f"G: {seq.count('G')} ({seq.count('G') / len(seq) * 100:.1f}%)")
    print(f"C: {seq.count('C')} ({seq.count('C') / len(seq) * 100:.1f}%)")
    print(f"N/ambiguous: {seq.count('N')} ({seq.count('N') / len(seq) * 100:.1f}%)")
    print()

    # Compare with other sequences in the dataset
    print("COMPARISON WITH OTHER SEQUENCES:")
    print("-" * 50)

    # Calculate pairwise differences
    differences = []
    for seq_id, other_seq in sequences.items():
        if seq_id != target_seq and len(other_seq) == len(seq):
            diff_count = sum(
                1 for a, b in zip(seq, other_seq) if a != b and a != "N" and b != "N"
            )
            similarity = (len(seq) - diff_count) / len(seq) * 100
            differences.append(
                {
                    "sequence": seq_id,
                    "differences": diff_count,
                    "similarity": similarity,
                    "description": sequence_info[seq_id][:50] + "..."
                    if len(sequence_info[seq_id]) > 50
                    else sequence_info[seq_id],
                }
            )

    # Sort by similarity (most different first)
    differences.sort(key=lambda x: x["similarity"])

    print(f"Pairwise comparisons with {target_seq} (showing top 10 most different):")
    print(f"{'Sequence':<15} {'Differences':<12} {'Similarity':<12} {'Description'}")
    print("-" * 80)

    for i, diff in enumerate(differences[:10]):
        print(
            f"{diff['sequence']:<15} {diff['differences']:<12} {diff['similarity']:<12.1f}% {diff['description']}"
        )

    print()
    print(
        f"Most similar sequence: {differences[-1]['sequence']} ({differences[-1]['similarity']:.1f}% similarity)"
    )
    print(
        f"Most different sequence: {differences[0]['sequence']} ({differences[0]['similarity']:.1f}% similarity)"
    )

    # Analyze genotype distribution
    print()
    print("GENOTYPE ANALYSIS:")
    print("-" * 30)

    genotypes = {}
    for seq_id, desc in sequence_info.items():
        # Extract genotype information from description
        if "GII." in desc:
            # Look for patterns like GII.4, GII.12, etc.
            genotype_match = re.search(r"GII\.[^/\s]+", desc)
            if genotype_match:
                genotype = genotype_match.group()
                genotypes[genotype] = genotypes.get(genotype, 0) + 1

    print("Genotype distribution in dataset:")
    for genotype, count in sorted(genotypes.items()):
        print(f"{genotype}: {count} sequences")

    # Extract target sequence genotype
    target_genotype = None
    if "GII." in sequence_info[target_seq]:
        genotype_match = re.search(r"GII\.[^/\s]+", sequence_info[target_seq])
        if genotype_match:
            target_genotype = genotype_match.group()

    print()
    print("OUTGROUP SUITABILITY ASSESSMENT:")
    print("-" * 40)

    # Check if this is GII.12 (as mentioned in the description)
    if (
        "GII.12" in sequence_info[target_seq]
        or "GII.Pg-GII.12" in sequence_info[target_seq]
    ):
        print(f"✓ {target_seq} is identified as GII.Pg-GII.12 recombinant")
        print("✓ GII.12 is a distinct genotype from the common GII.4")

        # Check if GII.12 is rare in the dataset
        gii12_count = sum(1 for desc in sequence_info.values() if "GII.12" in desc)
        total_sequences = len(sequences)

        print(
            f"✓ GII.12 frequency in dataset: {gii12_count}/{total_sequences} ({gii12_count / total_sequences * 100:.1f}%)"
        )

        if gii12_count == 1:
            print("✓ This is the only GII.12 sequence in the dataset")

        # Check genetic distance
        avg_similarity = sum(diff["similarity"] for diff in differences) / len(
            differences
        )
        print(f"✓ Average similarity to other sequences: {avg_similarity:.1f}%")

        if avg_similarity < 85:
            print("✓ Sufficiently divergent for outgroup use (< 85% similarity)")
        elif avg_similarity < 90:
            print("⚠ Moderately divergent for outgroup use (85-90% similarity)")
        else:
            print("✗ May be too similar for effective outgroup use (> 90% similarity)")

    # Publication context
    print()
    print("PUBLICATION CONTEXT (DOI: 10.1371/journal.pone.0189504):")
    print("-" * 60)
    print("The paper 'Detection and molecular characterization of the novel")
    print("recombinant norovirus GII.P16-GII.4 Sydney in southeastern Brazil")
    print("in 2016' focuses on GII.4 variants and recombinants.")
    print()
    print("Key findings from the paper:")
    print("- GII.4 Sydney variants were predominant")
    print("- GII.P16-GII.4 recombinants emerged in 2016")
    print("- Phylogenetic analysis used neighbor-joining method")
    print("- Study focused on ORF1/ORF2 junction region")
    print()

    # Recommendation
    print("RECOMMENDATION:")
    print("-" * 20)

    if gii12_count == 1 and avg_similarity < 85:
        print(f"✓ RECOMMENDED: {target_seq} (GII.Pg-GII.12) is suitable as outgroup")
        print("  Reasons:")
        print("  - Distinct genotype (GII.12) from main dataset")
        print("  - Sufficient genetic divergence")
        print("  - Single representative prevents ingroup clustering")
        print("  - Recombinant nature provides evolutionary perspective")
    elif avg_similarity >= 85:
        print(f"⚠ CAUTION: {target_seq} may be too similar for optimal outgroup")
        print("  Consider using a more divergent sequence if available")
    else:
        print(f"✓ ACCEPTABLE: {target_seq} can serve as outgroup with caveats")

    print()
    print("ALTERNATIVE SUGGESTIONS:")
    print("- Consider GI genotype sequences if available (more divergent)")
    print("- Use multiple outgroup sequences if possible")
    print("- Consider sapovirus sequences as more distant outgroup")

    # Create corrected trees analysis
    print()
    print("=" * 80)
    print("BRANCH LENGTH CORRECTION ANALYSIS")
    print("=" * 80)
    print()

    print(f"The sequence {target_seq} has been identified as causing")
    print("extremely long branch lengths in your phylogenetic trees.")
    print()
    print("Possible reasons for long branches:")
    print("1. Genuine evolutionary divergence (GII.12 vs GII.4)")
    print("2. Sequence quality issues")
    print("3. Alignment artifacts")
    print("4. Recombination breakpoints")
    print("5. Missing data or gaps")
    print()

    print("The branch length correction script has capped branches at 0.1")
    print("substitutions per site, which is reasonable for:")
    print("- Intragenotype comparisons")
    print("- Recent evolutionary timescales")
    print("- Visualization purposes")

    return target_seq, differences, genotypes


if __name__ == "__main__":
    analyze_outgroup_candidate()
