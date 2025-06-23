#!/usr/bin/env python3
"""
Additional Norovirus Sequence Analysis - Quality Assessment
"""

from Bio import SeqIO
import re
from collections import Counter
import pandas as pd


def check_sequence_quality(fasta_file):
    """Check sequence quality and potential issues"""

    sequences = []
    issues = []

    with open(fasta_file, "r") as f:
        for i, record in enumerate(SeqIO.parse(f, "fasta")):
            seq_str = str(record.seq).upper()

            # Basic sequence info
            seq_info = {
                "index": i + 1,
                "accession": record.id,
                "length": len(seq_str),
                "sequence": seq_str,
            }

            # Check for ambiguous nucleotides
            standard_nt = set("ATGC")
            seq_nt = set(seq_str)
            ambiguous_nt = seq_nt - standard_nt

            if ambiguous_nt:
                issues.append(
                    f"Sequence {i + 1} ({record.id}): Contains ambiguous nucleotides: {ambiguous_nt}"
                )
                seq_info["ambiguous_nucleotides"] = ambiguous_nt
            else:
                seq_info["ambiguous_nucleotides"] = set()

            # Check for unusual length (outliers)
            seq_info["is_outlier"] = False  # Will be determined later

            sequences.append(seq_info)

    # Determine outliers based on length
    lengths = [seq["length"] for seq in sequences]
    mean_length = sum(lengths) / len(lengths)
    std_length = (sum((x - mean_length) ** 2 for x in lengths) / len(lengths)) ** 0.5

    for seq in sequences:
        if abs(seq["length"] - mean_length) > 2 * std_length:
            seq["is_outlier"] = True
            issues.append(
                f"Sequence {seq['index']} ({seq['accession']}): "
                f"Length outlier ({seq['length']} bp, mean={mean_length:.1f})"
            )

    return sequences, issues


def analyze_start_end_sequences(sequences):
    """Analyze sequence start and end patterns"""

    start_sequences = []
    end_sequences = []

    for seq in sequences:
        seq_str = seq["sequence"]
        start_sequences.append(seq_str[:20])  # First 20 bp
        end_sequences.append(seq_str[-20:])  # Last 20 bp

    # Count common start/end patterns
    start_counts = Counter(start_sequences)
    end_counts = Counter(end_sequences)

    return start_counts, end_counts


def check_reading_frames(sequences):
    """Check for potential reading frame issues"""

    genetic_code = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    reading_frame_analysis = []

    for seq in sequences:
        seq_str = seq["sequence"]

        # Check all three reading frames
        for frame in range(3):
            codons = []
            stop_codons = 0

            for i in range(frame, len(seq_str) - 2, 3):
                codon = seq_str[i : i + 3]
                if len(codon) == 3:
                    codons.append(codon)
                    if genetic_code.get(codon, "X") == "*":
                        stop_codons += 1

            reading_frame_analysis.append(
                {
                    "accession": seq["accession"],
                    "frame": frame + 1,
                    "codons": len(codons),
                    "stop_codons": stop_codons,
                    "stop_percentage": (stop_codons / len(codons) * 100)
                    if codons
                    else 0,
                }
            )

    return reading_frame_analysis


def main():
    """Main analysis function"""
    fasta_file = (
        "/Users/berksakalli/Projects/automated-window-sliding/norovirus_sequences.fasta"
    )

    print("🔍 Additional Norovirus Sequence Quality Analysis")
    print("=" * 55)

    # Check sequence quality
    print("\n🧪 Checking sequence quality...")
    sequences, issues = check_sequence_quality(fasta_file)

    if issues:
        print(f"⚠️  Found {len(issues)} potential issues:")
        for issue in issues:
            print(f"   • {issue}")
    else:
        print("✅ No quality issues detected")

    # Analyze start/end patterns
    print(f"\n🔬 Analyzing sequence start/end patterns...")
    start_counts, end_counts = analyze_start_end_sequences(sequences)

    print(f"📍 Most common start sequences (first 20 bp):")
    for seq, count in start_counts.most_common(3):
        print(f"   • {seq}: {count} sequences ({count / len(sequences) * 100:.1f}%)")

    print(f"📍 Most common end sequences (last 20 bp):")
    for seq, count in end_counts.most_common(3):
        print(f"   • {seq}: {count} sequences ({count / len(sequences) * 100:.1f}%)")

    # Check reading frames
    print(f"\n🧬 Analyzing reading frames...")
    rf_analysis = check_reading_frames(sequences)

    # Create DataFrame for easier analysis
    df_rf = pd.DataFrame(rf_analysis)

    # Find the best reading frame for each sequence (least stop codons)
    best_frames = df_rf.groupby("accession")["stop_percentage"].idxmin()
    best_rf_data = df_rf.loc[best_frames]

    print(f"📊 Reading frame analysis summary:")
    print(
        f"   • Average stop codons in best frame: {best_rf_data['stop_percentage'].mean():.1f}%"
    )
    print(
        f"   • Range of stop codons: {best_rf_data['stop_percentage'].min():.1f}% - {best_rf_data['stop_percentage'].max():.1f}%"
    )

    # Identify sequences with high stop codon content
    high_stop_sequences = best_rf_data[best_rf_data["stop_percentage"] > 5]
    if not high_stop_sequences.empty:
        print(f"⚠️  Sequences with >5% stop codons in best frame:")
        for _, row in high_stop_sequences.iterrows():
            print(
                f"   • {row['accession']}: {row['stop_percentage']:.1f}% stop codons (frame {row['frame']})"
            )

    # Count sequences that start with ATG (start codon)
    atg_start_count = sum(1 for seq in sequences if seq["sequence"].startswith("ATG"))
    print(
        f"\n🎯 Sequences starting with ATG (start codon): {atg_start_count}/{len(sequences)} ({atg_start_count / len(sequences) * 100:.1f}%)"
    )

    # Save detailed analysis
    df_sequences = pd.DataFrame(
        [
            {
                "Accession": seq["accession"],
                "Length": seq["length"],
                "Is_Outlier": seq["is_outlier"],
                "Ambiguous_Nucleotides": len(seq["ambiguous_nucleotides"]),
                "Starts_with_ATG": seq["sequence"].startswith("ATG"),
                "First_20bp": seq["sequence"][:20],
                "Last_20bp": seq["sequence"][-20:],
            }
            for seq in sequences
        ]
    )

    quality_file = "/Users/berksakalli/Projects/automated-window-sliding/sequence_quality_analysis.csv"
    df_sequences.to_csv(quality_file, index=False)

    rf_file = "/Users/berksakalli/Projects/automated-window-sliding/reading_frame_analysis.csv"
    df_rf.to_csv(rf_file, index=False)

    print(f"\n💾 Detailed analysis saved to:")
    print(f"   • Quality analysis: {quality_file}")
    print(f"   • Reading frame analysis: {rf_file}")

    print(f"\n✅ Quality analysis complete!")

    # Summary recommendations
    print(f"\n💡 Recommendations for alignment:")

    outlier_count = sum(1 for seq in sequences if seq["is_outlier"])
    if outlier_count > 0:
        print(f"   • Consider reviewing {outlier_count} length outlier(s)")

    if len(start_counts) == 1:
        print(f"   • All sequences have identical start regions - good for alignment")
    else:
        print(f"   • {len(start_counts)} different start patterns - may need trimming")

    if best_rf_data["stop_percentage"].mean() < 2:
        print(f"   • Low stop codon content suggests good sequence quality")
    else:
        print(
            f"   • Higher stop codon content may indicate frameshifts or sequencing errors"
        )

    print(f"   • Consider using codon-aware alignment methods if translating")


if __name__ == "__main__":
    main()
