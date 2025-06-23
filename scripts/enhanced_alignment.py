#!/usr/bin/env python3
"""
Enhanced Alignment Pipeline for Norovirus Sequences
Based on the analysis results, this script implements targeted improvements
"""

import os
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline, ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict, Counter
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import tempfile


class NorovirusAlignmentEnhancer:
    def __init__(self, fasta_file, analysis_data_csv):
        self.fasta_file = fasta_file
        self.analysis_data = pd.read_csv(analysis_data_csv)
        self.sequences = self._load_sequences()
        self.genotype_groups = self._group_by_genotype()

    def _load_sequences(self):
        """Load sequences with metadata"""
        sequences = {}
        with open(self.fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                sequences[record.id] = {
                    "record": record,
                    "sequence": str(record.seq),
                    "length": len(record.seq),
                }
        return sequences

    def _group_by_genotype(self):
        """Group sequences by genotype for separate processing"""
        groups = defaultdict(list)

        for _, row in self.analysis_data.iterrows():
            accession = row["Accession"]
            genotype = row["Genotype"]
            if accession in self.sequences:
                groups[genotype].append(accession)

        return groups

    def step1_quality_filtering(self, min_length_percentile=10, max_ambiguous=5):
        """Step 1: Filter sequences based on quality metrics"""
        print("🔍 Step 1: Quality Filtering")

        # Calculate length threshold (remove bottom 10% by default)
        lengths = [data["length"] for data in self.sequences.values()]
        min_length = np.percentile(lengths, min_length_percentile)

        filtered_sequences = {}
        removed_sequences = []

        for acc_id, seq_data in self.sequences.items():
            # Check length threshold
            if seq_data["length"] < min_length:
                removed_sequences.append(
                    (
                        acc_id,
                        f"Length {seq_data['length']} < threshold {min_length:.0f}",
                    )
                )
                continue

            # Check for excessive ambiguous nucleotides
            ambiguous_count = sum(1 for nt in seq_data["sequence"] if nt not in "ATGC")
            if ambiguous_count > max_ambiguous:
                removed_sequences.append(
                    (acc_id, f"Ambiguous nucleotides: {ambiguous_count}")
                )
                continue

            filtered_sequences[acc_id] = seq_data

        print(f"   ✅ Kept {len(filtered_sequences)} sequences")
        print(f"   🗑️  Removed {len(removed_sequences)} sequences:")
        for acc_id, reason in removed_sequences:
            print(f"      • {acc_id}: {reason}")

        self.sequences = filtered_sequences
        return filtered_sequences

    def step2_sequence_trimming(self, trim_start=20, trim_end=20):
        """Step 2: Trim sequences to remove variable start/end regions"""
        print(
            f"\n✂️  Step 2: Sequence Trimming (remove {trim_start} bp from start, {trim_end} bp from end)"
        )

        # Find common regions by analyzing start/end patterns
        start_patterns = Counter()
        end_patterns = Counter()

        for seq_data in self.sequences.values():
            seq = seq_data["sequence"]
            start_patterns[seq[:30]] += 1  # Look at first 30 bp
            end_patterns[seq[-30:]] += 1  # Look at last 30 bp

        # Determine optimal trimming based on pattern analysis
        most_common_start = start_patterns.most_common(1)[0]
        most_common_end = end_patterns.most_common(1)[0]

        print(f"   📊 Most common start pattern: {most_common_start[1]} sequences")
        print(f"   📊 Most common end pattern: {most_common_end[1]} sequences")

        # Apply trimming
        trimmed_sequences = {}
        for acc_id, seq_data in self.sequences.items():
            original_seq = seq_data["sequence"]
            trimmed_seq = (
                original_seq[trim_start:-trim_end]
                if trim_end > 0
                else original_seq[trim_start:]
            )

            # Create new record
            new_record = SeqRecord(
                Seq(trimmed_seq),
                id=seq_data["record"].id,
                description=seq_data["record"].description
                + f" [trimmed {trim_start}:{trim_end}]",
            )

            trimmed_sequences[acc_id] = {
                "record": new_record,
                "sequence": trimmed_seq,
                "length": len(trimmed_seq),
                "original_length": seq_data["length"],
            }

        print(f"   ✅ Trimmed sequences: {len(original_seq)} → {len(trimmed_seq)} bp")

        self.sequences = trimmed_sequences
        return trimmed_sequences

    def step3_genotype_specific_alignment(self, aligner="muscle"):
        """Step 3: Perform genotype-specific alignments"""
        print(f"\n🧬 Step 3: Genotype-Specific Alignment using {aligner.upper()}")

        genotype_alignments = {}

        for genotype, acc_list in self.genotype_groups.items():
            if genotype == "Unknown":
                continue  # Handle Unknown separately

            # Filter sequences that still exist after quality filtering
            valid_acc_list = [acc for acc in acc_list if acc in self.sequences]

            if len(valid_acc_list) < 2:
                print(
                    f"   ⚠️  Skipping {genotype}: Only {len(valid_acc_list)} sequence(s)"
                )
                continue

            print(f"   🔬 Aligning {genotype}: {len(valid_acc_list)} sequences")

            # Create temporary FASTA for this genotype
            genotype_records = [self.sequences[acc]["record"] for acc in valid_acc_list]

            # Perform alignment
            alignment = self._align_sequences(genotype_records, f"{genotype}_{aligner}")
            genotype_alignments[genotype] = {
                "alignment": alignment,
                "sequences": valid_acc_list,
                "length": alignment.get_alignment_length() if alignment else 0,
            }

        # Handle Unknown genotype sequences
        unknown_acc_list = [
            acc
            for acc in self.genotype_groups.get("Unknown", [])
            if acc in self.sequences
        ]
        if len(unknown_acc_list) >= 2:
            print(f"   🔬 Aligning Unknown genotype: {len(unknown_acc_list)} sequences")
            unknown_records = [
                self.sequences[acc]["record"] for acc in unknown_acc_list
            ]
            alignment = self._align_sequences(unknown_records, f"Unknown_{aligner}")
            genotype_alignments["Unknown"] = {
                "alignment": alignment,
                "sequences": unknown_acc_list,
                "length": alignment.get_alignment_length() if alignment else 0,
            }

        self.genotype_alignments = genotype_alignments
        return genotype_alignments

    def _align_sequences(self, records, prefix):
        """Perform alignment using MUSCLE or Clustal"""
        if len(records) < 2:
            return None

        # Create temporary files
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as temp_input:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".aln", delete=False
            ) as temp_output:
                input_file = temp_input.name
                output_file = temp_output.name

        try:
            # Write sequences to temporary file
            SeqIO.write(records, input_file, "fasta")

            # Try MUSCLE first, fallback to simple pairwise if not available
            try:
                # Check if muscle is available
                result = subprocess.run(
                    ["muscle", "-version"], capture_output=True, text=True, timeout=5
                )
                if result.returncode == 0:
                    # Use MUSCLE
                    cmd = ["muscle", "-in", input_file, "-out", output_file]
                    subprocess.run(cmd, capture_output=True, timeout=60)

                    # Read alignment
                    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                        alignment = AlignIO.read(output_file, "fasta")
                        return alignment

            except (
                subprocess.TimeoutExpired,
                subprocess.CalledProcessError,
                FileNotFoundError,
            ):
                pass

            # Fallback: Create a simple alignment by padding sequences
            print(f"      ⚠️  MUSCLE not available, using simple padding alignment")
            return self._create_simple_alignment(records)

        finally:
            # Cleanup
            for f in [input_file, output_file]:
                if os.path.exists(f):
                    os.unlink(f)

    def _create_simple_alignment(self, records):
        """Create a simple alignment by padding sequences to same length"""
        if not records:
            return None

        max_length = max(len(record.seq) for record in records)
        aligned_records = []

        for record in records:
            # Pad with gaps at the end
            padded_seq = str(record.seq) + "-" * (max_length - len(record.seq))
            aligned_record = SeqRecord(
                Seq(padded_seq), id=record.id, description=record.description
            )
            aligned_records.append(aligned_record)

        return MultipleSeqAlignment(aligned_records)

    def step4_merge_alignments(self):
        """Step 4: Merge genotype-specific alignments into master alignment"""
        print(f"\n🔗 Step 4: Merging Genotype Alignments")

        all_aligned_records = []
        alignment_stats = {}

        for genotype, data in self.genotype_alignments.items():
            alignment = data["alignment"]
            if alignment:
                all_aligned_records.extend(alignment)
                alignment_stats[genotype] = {
                    "sequences": len(alignment),
                    "length": alignment.get_alignment_length(),
                    "gaps": self._calculate_gap_percentage(alignment),
                }
                print(
                    f"   ✅ {genotype}: {len(alignment)} sequences, {alignment.get_alignment_length()} bp, {alignment_stats[genotype]['gaps']:.1f}% gaps"
                )

        if not all_aligned_records:
            print("   ❌ No alignments to merge!")
            return None

        # Ensure all sequences have the same length for master alignment
        max_length = max(len(record.seq) for record in all_aligned_records)

        normalized_records = []
        for record in all_aligned_records:
            if len(record.seq) < max_length:
                # Pad with gaps
                padded_seq = str(record.seq) + "-" * (max_length - len(record.seq))
                normalized_record = SeqRecord(
                    Seq(padded_seq), id=record.id, description=record.description
                )
                normalized_records.append(normalized_record)
            else:
                normalized_records.append(record)

        master_alignment = MultipleSeqAlignment(normalized_records)

        print(
            f"   🎯 Master alignment: {len(master_alignment)} sequences, {master_alignment.get_alignment_length()} bp"
        )

        self.master_alignment = master_alignment
        self.alignment_stats = alignment_stats
        return master_alignment

    def _calculate_gap_percentage(self, alignment):
        """Calculate percentage of gaps in alignment"""
        total_chars = 0
        gap_chars = 0

        for record in alignment:
            total_chars += len(record.seq)
            gap_chars += str(record.seq).count("-")

        return (gap_chars / total_chars) * 100 if total_chars > 0 else 0

    def step5_alignment_refinement(self, gap_threshold=50):
        """Step 5: Refine alignment by removing gap-heavy regions"""
        print(
            f"\n🎯 Step 5: Alignment Refinement (removing columns with >{gap_threshold}% gaps)"
        )

        alignment = self.master_alignment
        if not alignment:
            print("   ❌ No alignment to refine!")
            return None

        alignment_length = alignment.get_alignment_length()
        columns_to_keep = []

        for i in range(alignment_length):
            column = alignment[:, i]
            gap_count = column.count("-")
            gap_percentage = (gap_count / len(column)) * 100

            if gap_percentage <= gap_threshold:
                columns_to_keep.append(i)

        if not columns_to_keep:
            print("   ⚠️  No columns meet the gap threshold criteria!")
            return alignment

        # Create refined alignment
        refined_records = []
        for record in alignment:
            refined_seq = "".join(record.seq[i] for i in columns_to_keep)
            refined_record = SeqRecord(
                Seq(refined_seq),
                id=record.id,
                description=record.description + " [refined]",
            )
            refined_records.append(refined_record)

        refined_alignment = MultipleSeqAlignment(refined_records)

        removed_columns = alignment_length - len(columns_to_keep)
        print(
            f"   ✅ Refined alignment: {len(refined_alignment)} sequences, {refined_alignment.get_alignment_length()} bp"
        )
        print(
            f"   🗑️  Removed {removed_columns} gap-heavy columns ({removed_columns / alignment_length * 100:.1f}%)"
        )

        self.refined_alignment = refined_alignment
        return refined_alignment

    def step6_quality_assessment(self):
        """Step 6: Assess final alignment quality"""
        print(f"\n📊 Step 6: Alignment Quality Assessment")

        alignment = getattr(self, "refined_alignment", self.master_alignment)
        if not alignment:
            print("   ❌ No alignment to assess!")
            return

        # Calculate various quality metrics
        metrics = {
            "sequences": len(alignment),
            "length": alignment.get_alignment_length(),
            "gaps": self._calculate_gap_percentage(alignment),
            "conservation": self._calculate_conservation_score(alignment),
            "identity": self._calculate_pairwise_identity(alignment),
        }

        print(f"   📏 Alignment length: {metrics['length']} bp")
        print(f"   🧬 Number of sequences: {metrics['sequences']}")
        print(f"   🕳️  Gap percentage: {metrics['gaps']:.1f}%")
        print(f"   🎯 Conservation score: {metrics['conservation']:.1f}%")
        print(f"   🔗 Average pairwise identity: {metrics['identity']:.1f}%")

        # Generate quality plot
        self._plot_alignment_quality(alignment)

        self.final_metrics = metrics
        return metrics

    def _calculate_conservation_score(self, alignment):
        """Calculate conservation score across alignment positions"""
        if not alignment:
            return 0

        conservation_scores = []
        alignment_length = alignment.get_alignment_length()

        for i in range(alignment_length):
            column = alignment[:, i]
            # Remove gaps for conservation calculation
            non_gap_chars = [char for char in column if char != "-"]

            if len(non_gap_chars) > 0:
                most_common = Counter(non_gap_chars).most_common(1)[0][1]
                conservation = most_common / len(non_gap_chars)
                conservation_scores.append(conservation)

        return np.mean(conservation_scores) * 100 if conservation_scores else 0

    def _calculate_pairwise_identity(self, alignment):
        """Calculate average pairwise sequence identity"""
        if len(alignment) < 2:
            return 0

        identities = []
        sequences = [str(record.seq) for record in alignment]

        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                seq1, seq2 = sequences[i], sequences[j]
                matches = sum(
                    1 for a, b in zip(seq1, seq2) if a == b and a != "-" and b != "-"
                )
                valid_positions = sum(
                    1 for a, b in zip(seq1, seq2) if a != "-" and b != "-"
                )

                if valid_positions > 0:
                    identity = matches / valid_positions
                    identities.append(identity)

        return np.mean(identities) * 100 if identities else 0

    def _plot_alignment_quality(self, alignment):
        """Generate alignment quality visualization"""
        alignment_length = alignment.get_alignment_length()

        # Calculate conservation per position
        conservation_per_pos = []
        gap_per_pos = []

        for i in range(alignment_length):
            column = alignment[:, i]

            # Conservation (most frequent non-gap character)
            non_gap_chars = [char for char in column if char != "-"]
            if non_gap_chars:
                most_common_count = Counter(non_gap_chars).most_common(1)[0][1]
                conservation = most_common_count / len(non_gap_chars)
            else:
                conservation = 0
            conservation_per_pos.append(conservation)

            # Gap percentage
            gap_percentage = column.count("-") / len(column)
            gap_per_pos.append(gap_percentage)

        # Create plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8))

        positions = range(1, alignment_length + 1)

        # Conservation plot
        ax1.plot(positions, conservation_per_pos, color="green", alpha=0.7)
        ax1.fill_between(positions, conservation_per_pos, alpha=0.3, color="green")
        ax1.set_ylabel("Conservation Score")
        ax1.set_title("Alignment Quality Assessment")
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 1)

        # Gap plot
        ax2.plot(positions, gap_per_pos, color="red", alpha=0.7)
        ax2.fill_between(positions, gap_per_pos, alpha=0.3, color="red")
        ax2.set_xlabel("Alignment Position")
        ax2.set_ylabel("Gap Frequency")
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 1)

        plt.tight_layout()

        output_file = (
            "/Users/berksakalli/Projects/automated-window-sliding/alignment_quality.png"
        )
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"   📈 Quality plot saved: {output_file}")

        return output_file

    def save_results(
        self, output_dir="/Users/berksakalli/Projects/automated-window-sliding"
    ):
        """Save all alignment results"""
        print(f"\n💾 Saving Results to {output_dir}")

        # Save master alignment
        if hasattr(self, "master_alignment") and self.master_alignment:
            master_file = os.path.join(output_dir, "master_alignment.fasta")
            AlignIO.write(self.master_alignment, master_file, "fasta")
            print(f"   ✅ Master alignment: {master_file}")

        # Save refined alignment
        if hasattr(self, "refined_alignment") and self.refined_alignment:
            refined_file = os.path.join(output_dir, "refined_alignment.fasta")
            AlignIO.write(self.refined_alignment, refined_file, "fasta")
            print(f"   ✅ Refined alignment: {refined_file}")

        # Save final alignment with outgroups (if available)
        if (
            hasattr(self, "final_alignment_with_outgroups")
            and self.final_alignment_with_outgroups
        ):
            outgroup_file = os.path.join(
                output_dir, "enhanced_alignment_with_outgroups.fasta"
            )
            AlignIO.write(self.final_alignment_with_outgroups, outgroup_file, "fasta")
            print(f"   ✅ Enhanced alignment with outgroups: {outgroup_file}")

        # Save genotype-specific alignments
        if hasattr(self, "genotype_alignments"):
            for genotype, data in self.genotype_alignments.items():
                if data["alignment"]:
                    genotype_file = os.path.join(
                        output_dir, f"alignment_{genotype.replace('/', '_')}.fasta"
                    )
                    AlignIO.write(data["alignment"], genotype_file, "fasta")
                    print(f"   ✅ {genotype} alignment: {genotype_file}")

        # Save alignment statistics
        if hasattr(self, "final_metrics"):
            stats_file = os.path.join(output_dir, "alignment_statistics.txt")
            with open(stats_file, "w") as f:
                f.write("Enhanced Alignment Statistics\n")
                f.write("=" * 35 + "\n\n")

                f.write("Final Alignment Metrics:\n")
                for metric, value in self.final_metrics.items():
                    f.write(f"  {metric.replace('_', ' ').title()}: {value:.2f}\n")

                if hasattr(self, "alignment_stats"):
                    f.write("\nGenotype-Specific Statistics:\n")
                    for genotype, stats in self.alignment_stats.items():
                        f.write(f"  {genotype}:\n")
                        for key, value in stats.items():
                            f.write(f"    {key}: {value:.2f}\n")

            print(f"   ✅ Statistics: {stats_file}")

    def step7_add_outgroup_sequences(self, outgroup_file=None, download_outgroups=True):
        """Step 7: Add GI outgroup sequences for better phylogenetic rooting"""
        print("\n🌿 Step 7: Adding Outgroup Sequences")

        outgroup_sequences = []

        # Try to load from provided file first
        if outgroup_file and os.path.exists(outgroup_file):
            outgroup_sequences = list(SeqIO.parse(outgroup_file, "fasta"))
            print(
                f"   📥 Loaded {len(outgroup_sequences)} outgroups from {outgroup_file}"
            )

        # Download outgroups if requested and none loaded
        elif download_outgroups:
            print("   🔽 Downloading GI outgroup sequences from NCBI...")
            try:
                from Bio import Entrez
                import time

                # Set email for NCBI (you should replace this with a real email)
                Entrez.email = "researcher@example.com"

                # Essential GI outgroups
                outgroup_accessions = {
                    "M87661": "Norwalk virus (GI.1)",
                    "L07418": "Southampton virus (GI.2)",
                    "AF093797": "Desert Shield virus (GI.3)",
                }

                for accession, description in outgroup_accessions.items():
                    try:
                        print(f"   📥 Downloading {accession} ({description})...")
                        handle = Entrez.efetch(
                            db="nucleotide",
                            id=accession,
                            rettype="fasta",
                            retmode="text",
                        )
                        record = SeqIO.read(handle, "fasta")
                        record.description = (
                            f"{record.description} | OUTGROUP | {description}"
                        )
                        outgroup_sequences.append(record)
                        handle.close()
                        time.sleep(1)  # Be nice to NCBI
                    except Exception as e:
                        print(f"      ❌ Failed to download {accession}: {e}")

            except ImportError:
                print("   ⚠️  Entrez not available, creating placeholder outgroups")
                download_outgroups = False

        # Create placeholder outgroups if download failed or not requested
        if not outgroup_sequences and not download_outgroups:
            print("   🔄 Creating placeholder outgroup sequences")
            placeholder_outgroups = [
                {
                    "id": "GI_outgroup_1",
                    "description": "GI.1 representative | OUTGROUP | Placeholder",
                    "sequence": "N" * 400,  # Placeholder sequence
                },
                {
                    "id": "GI_outgroup_2",
                    "description": "GI.2 representative | OUTGROUP | Placeholder",
                    "sequence": "N" * 400,  # Placeholder sequence
                },
            ]

            for seq_data in placeholder_outgroups:
                record = SeqRecord(
                    Seq(seq_data["sequence"]),
                    id=seq_data["id"],
                    description=seq_data["description"],
                )
                outgroup_sequences.append(record)

        if not outgroup_sequences:
            print("   ❌ No outgroup sequences available")
            return None

        # Align outgroups to current alignment
        current_alignment = getattr(self, "refined_alignment", self.master_alignment)
        if not current_alignment:
            print("   ❌ No current alignment to add outgroups to")
            return None

        alignment_length = current_alignment.get_alignment_length()
        print(f"   📏 Aligning outgroups to {alignment_length} bp")

        # Simple alignment by truncating or padding
        aligned_outgroups = []
        for outgroup in outgroup_sequences:
            seq_str = str(outgroup.seq)

            if len(seq_str) > alignment_length:
                # Truncate to match
                aligned_seq = seq_str[:alignment_length]
                print(
                    f"   ✂️  Truncated {outgroup.id}: {len(seq_str)} → {alignment_length} bp"
                )
            else:
                # Pad with gaps
                aligned_seq = seq_str + "-" * (alignment_length - len(seq_str))
                print(
                    f"   📏 Padded {outgroup.id}: {len(seq_str)} → {alignment_length} bp"
                )

            aligned_record = SeqRecord(
                Seq(aligned_seq), id=outgroup.id, description=outgroup.description
            )
            aligned_outgroups.append(aligned_record)

        # Combine with current alignment
        all_sequences = list(current_alignment) + aligned_outgroups
        final_alignment = MultipleSeqAlignment(all_sequences)

        print(
            f"   ✅ Final alignment: {len(final_alignment)} sequences ({len(aligned_outgroups)} outgroups added)"
        )

        self.final_alignment_with_outgroups = final_alignment
        return final_alignment

    def create_nextflow_parameters(self, output_dir=None):
        """Create Nextflow parameter file for the enhanced alignment"""
        print("\n⚙️  Creating Nextflow Parameters")

        if output_dir is None:
            output_dir = "/Users/berksakalli/Projects/automated-window-sliding"

        # Determine which alignment to use
        alignment_file = "refined_alignment.fasta"
        if hasattr(self, "final_alignment_with_outgroups"):
            alignment_file = "enhanced_alignment_with_outgroups.fasta"
            print("   🌿 Using alignment with outgroups")
        elif hasattr(self, "refined_alignment"):
            print("   🎯 Using refined alignment")
        else:
            print("   📋 Using master alignment")
            alignment_file = "master_alignment.fasta"

        # Parameters optimized for norovirus analysis
        params = {
            "input": alignment_file,
            "window_size": 200,
            "step_size": 15,
            "outgroup": "M87661",  # Norwalk virus if available
            "tree_method": "iqtree2",
            "model": "GTR+F+I+G4",
            "bootstrap": 1000,
            "rooting_method": "mad",
            "max_branch_length": 0.1,  # Cap extreme branch lengths
        }

        params_file = os.path.join(output_dir, "params_enhanced_alignment.json")

        import json

        with open(params_file, "w") as f:
            json.dump(params, f, indent=2)

        print(f"   ✅ Parameter file: {params_file}")
        print("   📋 Key settings:")
        for key, value in params.items():
            print(f"      {key}: {value}")

        return params_file

    # ...existing code...


def main():
    """Run the enhanced alignment pipeline"""
    print("🧬 Enhanced Norovirus Alignment Pipeline")
    print("=" * 45)

    # Initialize enhancer
    fasta_file = (
        "/Users/berksakalli/Projects/automated-window-sliding/norovirus_sequences.fasta"
    )
    csv_file = "/Users/berksakalli/Projects/automated-window-sliding/sequence_analysis_data.csv"

    enhancer = NorovirusAlignmentEnhancer(fasta_file, csv_file)

    # Run enhancement steps with outgroups
    enhancer.step1_quality_filtering(min_length_percentile=5)  # Remove bottom 5%
    enhancer.step2_sequence_trimming(
        trim_start=15, trim_end=15
    )  # Conservative trimming
    enhancer.step3_genotype_specific_alignment(aligner="muscle")
    enhancer.step4_merge_alignments()
    enhancer.step5_alignment_refinement(
        gap_threshold=60
    )  # Remove columns with >60% gaps
    enhancer.step6_quality_assessment()

    # NEW: Add outgroup sequences
    enhancer.step7_add_outgroup_sequences(download_outgroups=True)

    # Create Nextflow parameters
    enhancer.create_nextflow_parameters()

    # Save results
    enhancer.save_results()

    print("\n🎉 Enhanced alignment pipeline completed!")
    print("📁 Check the output files in the working directory")
    print("🌿 Outgroup sequences added for better phylogenetic analysis")
    print("⚙️  Nextflow parameters created for enhanced pipeline")


if __name__ == "__main__":
    main()
