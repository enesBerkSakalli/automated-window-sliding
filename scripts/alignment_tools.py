#!/usr/bin/env python3
"""
Alignment Validation and Alternative Methods
Provides additional tools for alignment enhancement
"""

import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import numpy as np
import subprocess
import os
import tempfile


def install_alignment_tools():
    """Check and provide installation instructions for alignment tools"""
    print("🔧 Checking Alignment Tools Availability")
    print("=" * 45)

    tools_status = {}

    # Check MUSCLE
    try:
        result = subprocess.run(
            ["muscle", "-version"], capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            tools_status["MUSCLE"] = True
            print("   ✅ MUSCLE: Available")
        else:
            tools_status["MUSCLE"] = False
    except (FileNotFoundError, subprocess.TimeoutExpired):
        tools_status["MUSCLE"] = False
        print("   ❌ MUSCLE: Not available")

    # Check MAFFT
    try:
        result = subprocess.run(
            ["mafft", "--version"], capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            tools_status["MAFFT"] = True
            print("   ✅ MAFFT: Available")
        else:
            tools_status["MAFFT"] = False
    except (FileNotFoundError, subprocess.TimeoutExpired):
        tools_status["MAFFT"] = False
        print("   ❌ MAFFT: Not available")

    # Check Clustal Omega
    try:
        result = subprocess.run(
            ["clustalo", "--version"], capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            tools_status["Clustal Omega"] = True
            print("   ✅ Clustal Omega: Available")
        else:
            tools_status["Clustal Omega"] = False
    except (FileNotFoundError, subprocess.TimeoutExpired):
        tools_status["Clustal Omega"] = False
        print("   ❌ Clustal Omega: Not available")

    # Provide installation instructions if tools are missing
    missing_tools = [tool for tool, available in tools_status.items() if not available]

    if missing_tools:
        print(f"\n📦 Installation Instructions for Missing Tools:")
        print("   To install missing alignment tools on macOS:")
        print("   ```bash")
        print("   # Install via Homebrew")
        print("   brew install muscle")
        print("   brew install mafft")
        print("   brew install clustal-omega")
        print("   ```")
        print("\n   Or via Conda:")
        print("   ```bash")
        print("   conda install -c bioconda muscle")
        print("   conda install -c bioconda mafft")
        print("   conda install -c bioconda clustal-omega")
        print("   ```")

    return tools_status


def create_codon_aware_alignment(fasta_file, output_file):
    """Create codon-aware alignment for protein-coding sequences"""
    print("🧬 Creating Codon-Aware Alignment")

    # Load sequences
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    # Translate sequences to check reading frames
    best_frame_sequences = []

    for record in sequences:
        seq_str = str(record.seq)

        # Try all three reading frames
        min_stops = float("inf")
        best_frame = 0

        for frame in range(3):
            frame_seq = seq_str[frame:]
            # Make length divisible by 3
            frame_seq = frame_seq[: len(frame_seq) - len(frame_seq) % 3]

            if len(frame_seq) >= 3:
                # Count stop codons
                stop_count = 0
                for i in range(0, len(frame_seq), 3):
                    codon = frame_seq[i : i + 3]
                    if codon.upper() in ["TAA", "TAG", "TGA"]:
                        stop_count += 1

                if stop_count < min_stops:
                    min_stops = stop_count
                    best_frame = frame

        # Use best frame for this sequence
        best_seq = seq_str[best_frame:]
        best_seq = best_seq[: len(best_seq) - len(best_seq) % 3]  # Make divisible by 3

        new_record = SeqRecord(
            Seq(best_seq),
            id=record.id,
            description=record.description + f" [frame_{best_frame + 1}]",
        )
        best_frame_sequences.append(new_record)

    # Save frame-adjusted sequences
    frame_file = output_file.replace(".fasta", "_frame_adjusted.fasta")
    SeqIO.write(best_frame_sequences, frame_file, "fasta")
    print(f"   ✅ Frame-adjusted sequences saved: {frame_file}")

    return frame_file


def run_multiple_aligners(fasta_file, output_dir):
    """Run multiple alignment tools and compare results"""
    print("🔄 Running Multiple Alignment Tools")

    results = {}

    # Try MAFFT (usually fastest and good quality)
    try:
        print("   🧬 Running MAFFT...")
        mafft_output = os.path.join(output_dir, "alignment_mafft.fasta")
        cmd = ["mafft", "--auto", fasta_file]

        with open(mafft_output, "w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, timeout=300)

        if result.returncode == 0 and os.path.exists(mafft_output):
            alignment = AlignIO.read(mafft_output, "fasta")
            results["MAFFT"] = {
                "file": mafft_output,
                "alignment": alignment,
                "length": alignment.get_alignment_length(),
                "sequences": len(alignment),
            }
            print(
                f"      ✅ MAFFT completed: {len(alignment)} seqs, {alignment.get_alignment_length()} bp"
            )

    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"      ❌ MAFFT failed: {e}")

    # Try MUSCLE
    try:
        print("   💪 Running MUSCLE...")
        muscle_output = os.path.join(output_dir, "alignment_muscle.fasta")
        cmd = ["muscle", "-in", fasta_file, "-out", muscle_output]

        result = subprocess.run(cmd, capture_output=True, timeout=300)

        if result.returncode == 0 and os.path.exists(muscle_output):
            alignment = AlignIO.read(muscle_output, "fasta")
            results["MUSCLE"] = {
                "file": muscle_output,
                "alignment": alignment,
                "length": alignment.get_alignment_length(),
                "sequences": len(alignment),
            }
            print(
                f"      ✅ MUSCLE completed: {len(alignment)} seqs, {alignment.get_alignment_length()} bp"
            )

    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"      ❌ MUSCLE failed: {e}")

    # Try Clustal Omega
    try:
        print("   🌟 Running Clustal Omega...")
        clustalo_output = os.path.join(output_dir, "alignment_clustalo.fasta")
        cmd = ["clustalo", "-i", fasta_file, "-o", clustalo_output, "--outfmt=fasta"]

        result = subprocess.run(cmd, capture_output=True, timeout=300)

        if result.returncode == 0 and os.path.exists(clustalo_output):
            alignment = AlignIO.read(clustalo_output, "fasta")
            results["Clustal Omega"] = {
                "file": clustalo_output,
                "alignment": alignment,
                "length": alignment.get_alignment_length(),
                "sequences": len(alignment),
            }
            print(
                f"      ✅ Clustal Omega completed: {len(alignment)} seqs, {alignment.get_alignment_length()} bp"
            )

    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"      ❌ Clustal Omega failed: {e}")

    return results


def compare_alignments(alignment_results):
    """Compare different alignment results"""
    print("\n📊 Comparing Alignment Results")

    if not alignment_results:
        print("   ❌ No alignments to compare")
        return

    comparison_data = []

    for method, data in alignment_results.items():
        alignment = data["alignment"]

        # Calculate metrics
        gap_percentage = calculate_gap_percentage(alignment)
        conservation_score = calculate_conservation_score(alignment)

        comparison_data.append(
            {
                "Method": method,
                "Sequences": len(alignment),
                "Length": alignment.get_alignment_length(),
                "Gap_%": gap_percentage,
                "Conservation_%": conservation_score,
            }
        )

    # Create comparison table
    df = pd.DataFrame(comparison_data)
    print(f"\n   Alignment Comparison:")
    print(df.to_string(index=False, float_format="%.1f"))

    # Find best alignment based on criteria
    # Prefer: low gap percentage, high conservation, reasonable length
    df["Score"] = (df["Conservation_%"] * 0.6) - (df["Gap_%"] * 0.4)
    best_method = df.loc[df["Score"].idxmax(), "Method"]

    print(f"\n   🏆 Recommended alignment: {best_method}")

    return df, best_method


def calculate_gap_percentage(alignment):
    """Calculate percentage of gaps in alignment"""
    total_chars = sum(len(record.seq) for record in alignment)
    gap_chars = sum(str(record.seq).count("-") for record in alignment)
    return (gap_chars / total_chars) * 100 if total_chars > 0 else 0


def calculate_conservation_score(alignment):
    """Calculate conservation score"""
    if not alignment:
        return 0

    conservation_scores = []
    alignment_length = alignment.get_alignment_length()

    for i in range(alignment_length):
        column = alignment[:, i]
        non_gap_chars = [char for char in column if char != "-"]

        if len(non_gap_chars) > 1:
            from collections import Counter

            most_common_count = Counter(non_gap_chars).most_common(1)[0][1]
            conservation = most_common_count / len(non_gap_chars)
            conservation_scores.append(conservation)

    return np.mean(conservation_scores) * 100 if conservation_scores else 0


def create_alignment_report(fasta_file, analysis_csv, output_dir):
    """Create comprehensive alignment enhancement report"""
    print("📝 Creating Alignment Enhancement Report")

    # Load analysis data
    df = pd.read_csv(analysis_csv)

    report_file = os.path.join(output_dir, "alignment_enhancement_report.md")

    with open(report_file, "w") as f:
        f.write("# Norovirus Sequence Alignment Enhancement Report\n\n")

        f.write("## Dataset Overview\n")
        f.write(f"- **Total sequences:** {len(df)}\n")
        f.write(
            f"- **Length range:** {df['Length_bp'].min()} - {df['Length_bp'].max()} bp\n"
        )
        f.write(
            f"- **Mean length:** {df['Length_bp'].mean():.1f} ± {df['Length_bp'].std():.1f} bp\n"
        )
        f.write(
            f"- **Genotypes:** {df['Genotype'].nunique()} unique ({df['Genotype'].value_counts().head(3).to_dict()})\n\n"
        )

        f.write("## Identified Challenges\n")
        f.write(
            "Based on our sequence analysis, the following challenges were identified:\n\n"
        )

        # Length outliers
        outliers = df[
            (df["Length_bp"] < df["Length_bp"].quantile(0.1))
            | (df["Length_bp"] > df["Length_bp"].quantile(0.9))
        ]
        f.write(
            f"1. **Length outliers:** {len(outliers)} sequences with extreme lengths\n"
        )
        for _, row in outliers.iterrows():
            f.write(f"   - {row['Accession']}: {row['Length_bp']} bp\n")
        f.write("\n")

        # Genotype diversity
        f.write(
            f"2. **Genotype diversity:** {df['Genotype'].nunique()} different genotypes\n"
        )
        f.write(
            "   - This suggests high sequence divergence requiring careful alignment\n\n"
        )

        # GC content variation
        gc_range = df["GC_Content_%"].max() - df["GC_Content_%"].min()
        f.write(f"3. **GC content variation:** {gc_range:.1f}% range\n")
        f.write("   - May indicate compositional bias affecting alignment\n\n")

        f.write("## Enhancement Strategy\n")
        f.write("Our multi-step enhancement approach addresses these challenges:\n\n")

        f.write("### Step 1: Quality Filtering\n")
        f.write("- Remove sequences with extreme lengths (bottom 5%)\n")
        f.write("- Filter sequences with excessive ambiguous nucleotides\n")
        f.write("- Retain high-quality sequences for alignment\n\n")

        f.write("### Step 2: Sequence Trimming\n")
        f.write("- Remove variable start/end regions (15 bp each end)\n")
        f.write("- Focus on conserved core regions\n")
        f.write("- Reduce alignment artifacts from partial sequences\n\n")

        f.write("### Step 3: Genotype-Specific Alignment\n")
        f.write("- Align sequences within genotype groups first\n")
        f.write("- Reduces impact of high divergence between genotypes\n")
        f.write("- Improves local alignment accuracy\n\n")

        f.write("### Step 4: Alignment Merging\n")
        f.write("- Combine genotype-specific alignments\n")
        f.write("- Maintain within-group alignment quality\n")
        f.write("- Create comprehensive master alignment\n\n")

        f.write("### Step 5: Gap Refinement\n")
        f.write("- Remove columns with >60% gaps\n")
        f.write("- Focus on informative alignment regions\n")
        f.write("- Improve signal-to-noise ratio\n\n")

        f.write("### Step 6: Quality Assessment\n")
        f.write("- Calculate conservation scores\n")
        f.write("- Assess pairwise sequence identity\n")
        f.write("- Generate quality visualizations\n\n")

        f.write("## Recommended Tools\n")
        f.write("For optimal results, we recommend:\n\n")
        f.write("1. **MAFFT** - Fast and accurate for large datasets\n")
        f.write("2. **MUSCLE** - Good balance of speed and accuracy\n")
        f.write("3. **Clustal Omega** - High accuracy for divergent sequences\n\n")

        f.write("## Expected Outcomes\n")
        f.write("This enhancement approach should result in:\n\n")
        f.write("- **Improved conservation scores** in functional regions\n")
        f.write("- **Reduced gap artifacts** from sequence length variation\n")
        f.write("- **Better phylogenetic signal** for downstream analysis\n")
        f.write("- **Optimal window selection** for sliding window analysis\n\n")

        f.write("## Next Steps\n")
        f.write("After enhanced alignment:\n\n")
        f.write("1. **Validate alignment quality** using multiple metrics\n")
        f.write("2. **Identify conserved regions** for primer design\n")
        f.write("3. **Set sliding window parameters** based on alignment length\n")
        f.write("4. **Consider phylogenetic analysis** to understand relationships\n")

    print(f"   ✅ Report saved: {report_file}")
    return report_file


def main():
    """Run alignment validation and enhancement tools"""
    print("🔧 Alignment Enhancement Tools & Validation")
    print("=" * 50)

    # File paths
    fasta_file = (
        "/Users/berksakalli/Projects/automated-window-sliding/norovirus_sequences.fasta"
    )
    analysis_csv = "/Users/berksakalli/Projects/automated-window-sliding/sequence_analysis_data.csv"
    output_dir = "/Users/berksakalli/Projects/automated-window-sliding"

    # Check tool availability
    tools_status = install_alignment_tools()

    # Create enhancement report
    create_alignment_report(fasta_file, analysis_csv, output_dir)

    # Create codon-aware version if requested
    print(f"\n🧬 Creating codon-aware alignment preparation...")
    codon_file = create_codon_aware_alignment(
        fasta_file, os.path.join(output_dir, "codon_aware_sequences.fasta")
    )

    # Run multiple aligners if available
    available_tools = [tool for tool, status in tools_status.items() if status]
    if available_tools:
        print(f"\n🔄 Running available alignment tools: {', '.join(available_tools)}")
        alignment_results = run_multiple_aligners(fasta_file, output_dir)

        if alignment_results:
            comparison_df, best_method = compare_alignments(alignment_results)

            # Save comparison results
            comparison_file = os.path.join(output_dir, "alignment_comparison.csv")
            comparison_df.to_csv(comparison_file, index=False)
            print(f"   📊 Comparison saved: {comparison_file}")
    else:
        print(
            f"\n⚠️  No alignment tools available. Please install MAFFT, MUSCLE, or Clustal Omega."
        )

    print(f"\n✅ Alignment enhancement tools completed!")


if __name__ == "__main__":
    main()
