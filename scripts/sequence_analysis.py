#!/usr/bin/env python3
"""
Norovirus Sequence Analysis Script
Performs comprehensive analysis of FASTA sequences before alignment
"""

import re
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def parse_fasta_info(fasta_file):
    """Parse FASTA file and extract sequence information"""
    sequences = []

    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            # Extract metadata from header
            header_parts = record.description.split()
            accession = header_parts[0] if header_parts else "Unknown"

            # Try to extract genotype information
            genotype = "Unknown"
            if "GII" in record.description:
                genotype_match = re.search(r"GII\.[P\d]+-GII\.\d+", record.description)
                if genotype_match:
                    genotype = genotype_match.group()

            sequences.append(
                {
                    "accession": accession,
                    "description": record.description,
                    "genotype": genotype,
                    "sequence": str(record.seq),
                    "length": len(record.seq),
                    "gc_content": gc_fraction(record.seq) * 100,
                    "molecular_weight": molecular_weight(record.seq, seq_type="DNA"),
                }
            )

    return sequences


def analyze_nucleotide_composition(sequences):
    """Analyze nucleotide composition across all sequences"""
    all_nucleotides = "".join([seq["sequence"] for seq in sequences])
    nt_counts = Counter(all_nucleotides.upper())

    total_nt = sum(nt_counts.values())
    nt_percentages = {nt: (count / total_nt) * 100 for nt, count in nt_counts.items()}

    return nt_counts, nt_percentages


def analyze_sequence_lengths(sequences):
    """Analyze sequence length distribution"""
    lengths = [seq["length"] for seq in sequences]

    return {
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": np.mean(lengths),
        "median_length": np.median(lengths),
        "std_length": np.std(lengths),
        "lengths": lengths,
    }


def analyze_gc_content(sequences):
    """Analyze GC content distribution"""
    gc_contents = [seq["gc_content"] for seq in sequences]

    return {
        "min_gc": min(gc_contents),
        "max_gc": max(gc_contents),
        "mean_gc": np.mean(gc_contents),
        "median_gc": np.median(gc_contents),
        "std_gc": np.std(gc_contents),
        "gc_contents": gc_contents,
    }


def analyze_genotypes(sequences):
    """Analyze genotype distribution"""
    genotypes = [seq["genotype"] for seq in sequences]
    genotype_counts = Counter(genotypes)

    return genotype_counts


def find_conserved_regions(sequences, min_length=10):
    """Find conserved regions across sequences"""
    if not sequences:
        return []

    # Get the shortest sequence length for comparison
    min_seq_length = min(len(seq["sequence"]) for seq in sequences)

    conserved_regions = []

    # Check for conserved regions of varying lengths
    for region_length in range(min_length, min(50, min_seq_length + 1)):
        for start_pos in range(min_seq_length - region_length + 1):
            # Extract the region from the first sequence
            reference_region = sequences[0]["sequence"][
                start_pos : start_pos + region_length
            ]

            # Check if this region is conserved across all sequences
            is_conserved = True
            for seq in sequences[1:]:
                if len(seq["sequence"]) >= start_pos + region_length:
                    seq_region = seq["sequence"][start_pos : start_pos + region_length]
                    if seq_region != reference_region:
                        is_conserved = False
                        break
                else:
                    is_conserved = False
                    break

            if is_conserved:
                conserved_regions.append(
                    {
                        "start": start_pos,
                        "end": start_pos + region_length,
                        "length": region_length,
                        "sequence": reference_region,
                    }
                )

    # Remove overlapping regions, keeping the longest ones
    conserved_regions.sort(key=lambda x: x["length"], reverse=True)
    filtered_regions = []

    for region in conserved_regions:
        overlap = False
        for existing in filtered_regions:
            if region["start"] < existing["end"] and region["end"] > existing["start"]:
                overlap = True
                break
        if not overlap:
            filtered_regions.append(region)

    return filtered_regions


def generate_plots(sequences, length_stats, gc_stats, nt_percentages, genotype_counts):
    """Generate visualization plots"""

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle("Norovirus Sequence Analysis", fontsize=16, fontweight="bold")

    # 1. Sequence length distribution
    axes[0, 0].hist(
        length_stats["lengths"], bins=20, alpha=0.7, color="skyblue", edgecolor="black"
    )
    axes[0, 0].axvline(
        length_stats["mean_length"],
        color="red",
        linestyle="--",
        label=f"Mean: {length_stats['mean_length']:.1f}",
    )
    axes[0, 0].set_xlabel("Sequence Length (bp)")
    axes[0, 0].set_ylabel("Frequency")
    axes[0, 0].set_title("Sequence Length Distribution")
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # 2. GC content distribution
    axes[0, 1].hist(
        gc_stats["gc_contents"],
        bins=20,
        alpha=0.7,
        color="lightgreen",
        edgecolor="black",
    )
    axes[0, 1].axvline(
        gc_stats["mean_gc"],
        color="red",
        linestyle="--",
        label=f"Mean: {gc_stats['mean_gc']:.1f}%",
    )
    axes[0, 1].set_xlabel("GC Content (%)")
    axes[0, 1].set_ylabel("Frequency")
    axes[0, 1].set_title("GC Content Distribution")
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # 3. Nucleotide composition
    nucleotides = ["A", "T", "G", "C"]
    percentages = [nt_percentages.get(nt, 0) for nt in nucleotides]
    colors = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"]

    axes[0, 2].pie(
        percentages, labels=nucleotides, autopct="%1.1f%%", colors=colors, startangle=90
    )
    axes[0, 2].set_title("Overall Nucleotide Composition")

    # 4. Genotype distribution
    if len(genotype_counts) > 1:
        genotypes = list(genotype_counts.keys())
        counts = list(genotype_counts.values())

        axes[1, 0].bar(range(len(genotypes)), counts, color="coral", alpha=0.7)
        axes[1, 0].set_xlabel("Genotype")
        axes[1, 0].set_ylabel("Count")
        axes[1, 0].set_title("Genotype Distribution")
        axes[1, 0].set_xticks(range(len(genotypes)))
        axes[1, 0].set_xticklabels(genotypes, rotation=45, ha="right")
        axes[1, 0].grid(True, alpha=0.3)
    else:
        axes[1, 0].text(
            0.5,
            0.5,
            "Single genotype detected",
            ha="center",
            va="center",
            transform=axes[1, 0].transAxes,
            fontsize=12,
        )
        axes[1, 0].set_title("Genotype Distribution")

    # 5. Sequence length vs GC content scatter plot
    lengths = [seq["length"] for seq in sequences]
    gc_contents = [seq["gc_content"] for seq in sequences]

    axes[1, 1].scatter(lengths, gc_contents, alpha=0.7, color="purple", s=50)
    axes[1, 1].set_xlabel("Sequence Length (bp)")
    axes[1, 1].set_ylabel("GC Content (%)")
    axes[1, 1].set_title("Sequence Length vs GC Content")
    axes[1, 1].grid(True, alpha=0.3)

    # 6. Summary statistics table
    axes[1, 2].axis("off")

    summary_data = [
        ["Total Sequences", len(sequences)],
        ["Min Length", f"{length_stats['min_length']} bp"],
        ["Max Length", f"{length_stats['max_length']} bp"],
        ["Mean Length", f"{length_stats['mean_length']:.1f} bp"],
        ["Mean GC Content", f"{gc_stats['mean_gc']:.1f}%"],
        ["Unique Genotypes", len(genotype_counts)],
    ]

    table = axes[1, 2].table(
        cellText=summary_data,
        colLabels=["Metric", "Value"],
        cellLoc="left",
        loc="center",
        colWidths=[0.6, 0.4],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style the table
    for i in range(len(summary_data) + 1):
        for j in range(2):
            cell = table[(i, j)]
            if i == 0:  # Header row
                cell.set_facecolor("#4472C4")
                cell.set_text_props(weight="bold", color="white")
            else:
                cell.set_facecolor("#F2F2F2" if i % 2 == 0 else "white")

    axes[1, 2].set_title("Summary Statistics")

    plt.tight_layout()
    return fig


def main():
    """Main analysis function"""
    fasta_file = (
        "/Users/berksakalli/Projects/automated-window-sliding/norovirus_sequences.fasta"
    )

    print("🧬 Norovirus Sequence Analysis Report")
    print("=" * 50)

    # Parse sequences
    print("\n📁 Loading sequences...")
    sequences = parse_fasta_info(fasta_file)
    print(f"✅ Loaded {len(sequences)} sequences")

    # Basic statistics
    print("\n📊 Basic Statistics:")
    print(f"   • Total sequences: {len(sequences)}")

    # Analyze sequence lengths
    length_stats = analyze_sequence_lengths(sequences)
    print(f"\n📏 Sequence Length Analysis:")
    print(f"   • Min length: {length_stats['min_length']} bp")
    print(f"   • Max length: {length_stats['max_length']} bp")
    print(f"   • Mean length: {length_stats['mean_length']:.1f} bp")
    print(f"   • Median length: {length_stats['median_length']:.1f} bp")
    print(f"   • Standard deviation: {length_stats['std_length']:.1f} bp")

    # Analyze GC content
    gc_stats = analyze_gc_content(sequences)
    print(f"\n🧪 GC Content Analysis:")
    print(f"   • Min GC content: {gc_stats['min_gc']:.1f}%")
    print(f"   • Max GC content: {gc_stats['max_gc']:.1f}%")
    print(f"   • Mean GC content: {gc_stats['mean_gc']:.1f}%")
    print(f"   • Median GC content: {gc_stats['median_gc']:.1f}%")
    print(f"   • Standard deviation: {gc_stats['std_gc']:.1f}%")

    # Analyze nucleotide composition
    nt_counts, nt_percentages = analyze_nucleotide_composition(sequences)
    print(f"\n🔬 Nucleotide Composition:")
    for nt in ["A", "T", "G", "C"]:
        if nt in nt_percentages:
            print(
                f"   • {nt}: {nt_percentages[nt]:.1f}% ({nt_counts[nt]:,} nucleotides)"
            )

    # Analyze genotypes
    genotype_counts = analyze_genotypes(sequences)
    print(f"\n🦠 Genotype Distribution:")
    for genotype, count in genotype_counts.most_common():
        print(
            f"   • {genotype}: {count} sequences ({count / len(sequences) * 100:.1f}%)"
        )

    # Find conserved regions
    print(f"\n🔍 Searching for conserved regions...")
    conserved_regions = find_conserved_regions(sequences, min_length=15)

    if conserved_regions:
        print(f"   • Found {len(conserved_regions)} conserved regions:")
        for i, region in enumerate(conserved_regions[:5], 1):  # Show top 5
            print(
                f"     {i}. Position {region['start']}-{region['end']} "
                f"({region['length']} bp): {region['sequence']}"
            )
        if len(conserved_regions) > 5:
            print(f"     ... and {len(conserved_regions) - 5} more regions")
    else:
        print("   • No highly conserved regions found with minimum length of 15 bp")

    # Generate and save plots
    print(f"\n📈 Generating visualizations...")
    fig = generate_plots(
        sequences, length_stats, gc_stats, nt_percentages, genotype_counts
    )

    output_file = "/Users/berksakalli/Projects/automated-window-sliding/sequence_analysis_plots.png"
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"   • Plots saved to: {output_file}")

    # Save detailed data to CSV
    df_data = []
    for seq in sequences:
        df_data.append(
            {
                "Accession": seq["accession"],
                "Genotype": seq["genotype"],
                "Length_bp": seq["length"],
                "GC_Content_%": seq["gc_content"],
                "Molecular_Weight": seq["molecular_weight"],
                "Description": seq["description"],
            }
        )

    df = pd.DataFrame(df_data)
    csv_file = "/Users/berksakalli/Projects/automated-window-sliding/sequence_analysis_data.csv"
    df.to_csv(csv_file, index=False)
    print(f"   • Detailed data saved to: {csv_file}")

    print(f"\n✅ Analysis complete!")
    print(f"\n💡 Key Insights:")

    # Generate insights
    length_range = length_stats["max_length"] - length_stats["min_length"]
    if length_range > 100:
        print(f"   • High length variability detected ({length_range} bp range)")
        print(f"     This may affect alignment quality and should be considered")
    else:
        print(
            f"   • Sequences have relatively consistent lengths ({length_range} bp range)"
        )

    gc_range = gc_stats["max_gc"] - gc_stats["min_gc"]
    if gc_range > 10:
        print(f"   • Significant GC content variation ({gc_range:.1f}% range)")
        print(f"     This suggests potential sequence diversity")
    else:
        print(f"   • GC content is relatively consistent ({gc_range:.1f}% range)")

    if len(genotype_counts) > 1:
        print(f"   • Multiple genotypes detected ({len(genotype_counts)} unique)")
        print(f"     Consider genotype-specific analysis")

    if conserved_regions:
        avg_conserved_length = np.mean([r["length"] for r in conserved_regions])
        print(
            f"   • Conserved regions identified (avg length: {avg_conserved_length:.1f} bp)"
        )
        print(f"     These may be good targets for primer design")


if __name__ == "__main__":
    main()
