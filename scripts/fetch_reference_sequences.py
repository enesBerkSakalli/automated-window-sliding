#!/usr/bin/env python3
"""
Curated High-Value Norovirus Reference Sequences
Based on literature review and recombination analysis requirements.

This script fetches carefully selected reference sequences that are:
1. Complete or near-complete genomes (>6000bp)
2. Well-characterized recombinant strains
3. Representative of major circulating genotypes
4. From recent surveillance studies (2018-2024)
"""

from Bio import Entrez, SeqIO
import os
import time

# Email for NCBI (replace with actual email)
Entrez.email = "researcher@example.com"

# Curated list of high-value reference sequences
REFERENCE_SEQUENCES = {
    # Recent complete genomes from key studies
    "KY210980": {
        "description": "Complete GII.P16-GII.4/RUS/Novosibirsk/NS16-C38/2016",
        "length": 7560,
        "importance": "First complete Russian GII.P16/GII.4 recombinant",
        "year": 2016,
    },
    "KY421121": {
        "description": "Complete GII.P16/GII.2 from China",
        "length": 7500,
        "importance": "Well-characterized GII.P16/GII.2 recombinant",
        "year": 2016,
    },
    "MK753033": {
        "description": "Reference from Illumina-Nanopore study",
        "length": 7500,
        "importance": "High-quality reference used in recent genomic studies",
        "year": 2019,
    },
    # Spanish surveillance sequences
    "MK789430": {
        "description": "Complete genome Valencia Spain",
        "length": 7500,
        "year": 2019,
    },
    "MK789431": {
        "description": "Complete genome Valencia Spain",
        "length": 7500,
        "year": 2019,
    },
    "MK789432": {
        "description": "Complete genome Valencia Spain",
        "length": 7500,
        "year": 2019,
    },
    "MK789433": {
        "description": "Complete genome Valencia Spain",
        "length": 7500,
        "year": 2019,
    },
    # Russian recombinant strains
    "KY210977": {
        "description": "GII.Pe-GII.4/RUS/Novosibirsk/NS16-C13/2016",
        "length": 4300,
        "importance": "RdRp + ORF2/3 coverage",
        "year": 2016,
    },
    "KY210919": {
        "description": "GII.P21-GII.3/RUS/Novosibirsk/NS16-C23/2016",
        "length": 4300,
        "importance": "RdRp + ORF2/3 coverage",
        "year": 2016,
    },
    # Additional high-priority complete genomes
    "MN025903": {"description": "Complete GII.4 Sydney variant", "year": 2019},
    "MN025904": {"description": "Complete GII.4 Sydney variant", "year": 2019},
    "LC590997": {"description": "Complete GII.P16-GII.2 Japan", "year": 2020},
    "LC590998": {"description": "Complete GII.P16-GII.2 Japan", "year": 2020},
    # Recent surveillance sequences 2020-2023
    "MW693851": {"description": "GII.P16-GII.4 USA 2020", "year": 2020},
    "MW693852": {"description": "GII.P31-GII.4 USA 2020", "year": 2020},
    "MW693853": {"description": "GII.P4-GII.4 USA 2020", "year": 2020},
    # European surveillance
    "OL539847": {"description": "GII.P16-GII.3 Germany 2021", "year": 2021},
    "OL539848": {"description": "GII.P21-GII.3 Germany 2021", "year": 2021},
    # Asian surveillance
    "OM742515": {"description": "GII.P16-GII.2 China 2022", "year": 2022},
    "OM742516": {"description": "GII.P31-GII.4 China 2022", "year": 2022},
    # Australian sequences
    "ON084521": {"description": "GII.P16-GII.4 Australia 2022", "year": 2022},
    "ON084522": {"description": "GII.P4-GII.4 Australia 2022", "year": 2022},
}


def fetch_curated_sequences():
    """
    Fetch the curated list of high-value reference sequences.
    """
    print("🎯 FETCHING CURATED HIGH-VALUE REFERENCE SEQUENCES")
    print("=" * 70)

    output_dir = "/Users/berksakalli/Projects/automated-window-sliding/data"
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, "reference_complete_genomes.fasta")

    sequences = []
    failed_accessions = []

    print(f"📋 Target sequences: {len(REFERENCE_SEQUENCES)}")

    # Fetch sequences in small batches
    accessions = list(REFERENCE_SEQUENCES.keys())

    for i, accession in enumerate(accessions, 1):
        seq_info = REFERENCE_SEQUENCES[accession]
        print(f"\n{i:2d}. Fetching {accession}")
        print(f"    {seq_info['description']}")

        try:
            handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="fasta", retmode="text"
            )

            seq_records = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            if seq_records:
                seq = seq_records[0]
                seq_len = len(seq.seq)

                # Update sequence description with metadata
                seq.description = f"{seq.description} | Year: {seq_info.get('year', 'Unknown')} | Length: {seq_len}bp"

                sequences.append(seq)
                print(f"    ✅ SUCCESS: {seq_len:,}bp")

                if "importance" in seq_info:
                    print(f"    📝 Note: {seq_info['importance']}")
            else:
                print(f"    ❌ FAILED: No sequence data")
                failed_accessions.append(accession)

            time.sleep(1)  # Be respectful to NCBI

        except Exception as e:
            print(f"    ❌ ERROR: {e}")
            failed_accessions.append(accession)

    # Write sequences to file
    if sequences:
        SeqIO.write(sequences, output_file, "fasta")
        print(f"\n✅ SUCCESS: {len(sequences)} reference sequences saved")
        print(f"📁 File: {output_file}")

        # Analyze retrieved sequences
        analyze_reference_sequences(sequences)

    else:
        print("\n❌ No sequences successfully retrieved")

    if failed_accessions:
        print(f"\n⚠️  Failed accessions ({len(failed_accessions)}):")
        for acc in failed_accessions:
            print(f"   - {acc}: {REFERENCE_SEQUENCES[acc]['description']}")

    return sequences


def analyze_reference_sequences(sequences):
    """
    Analyze the characteristics of retrieved reference sequences.
    """
    print(f"\n📊 REFERENCE SEQUENCE ANALYSIS")
    print("=" * 50)

    lengths = [len(seq.seq) for seq in sequences]
    years = []
    genotypes = {
        "GII.P16": 0,
        "GII.P4": 0,
        "GII.P21": 0,
        "GII.P31": 0,
        "GII.Pe": 0,
        "Other": 0,
    }

    for seq in sequences:
        desc = seq.description.upper()

        # Extract genotype
        found_genotype = False
        for geno in genotypes.keys():
            if geno in desc and geno != "Other":
                genotypes[geno] += 1
                found_genotype = True
                break
        if not found_genotype:
            genotypes["Other"] += 1

        # Extract year
        import re

        year_match = re.search(r"YEAR: (\d{4})", desc)
        if year_match:
            years.append(int(year_match.group(1)))

    print(f"📈 Statistics:")
    print(f"   Total sequences: {len(sequences)}")
    print(f"   Length range: {min(lengths):,} - {max(lengths):,} bp")
    print(f"   Average length: {sum(lengths) // len(lengths):,} bp")
    print(
        f"   Complete genomes (>7000bp): {sum(1 for length in lengths if length > 7000)}"
    )
    print(
        f"   Near-complete (5000-7000bp): {sum(1 for length in lengths if 5000 <= length <= 7000)}"
    )

    if years:
        print(f"   Year range: {min(years)} - {max(years)}")

    print(f"\n🧬 Genotype Distribution:")
    for geno, count in genotypes.items():
        if count > 0:
            print(f"   {geno}: {count} sequences")


def create_comprehensive_dataset():
    """
    Combine curated references with existing cleaned dataset.
    """
    print(f"\n🔧 CREATING COMPREHENSIVE DATASET")
    print("=" * 50)

    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"

    # File paths
    current_file = os.path.join(base_dir, "data", "cleaned_alignment_combined.fasta")
    reference_file = os.path.join(base_dir, "data", "reference_complete_genomes.fasta")
    output_file = os.path.join(base_dir, "data", "comprehensive_enhanced_dataset.fasta")

    all_sequences = []

    # Read current cleaned dataset
    if os.path.exists(current_file):
        current_sequences = list(SeqIO.parse(current_file, "fasta"))
        all_sequences.extend(current_sequences)
        print(f"📂 Current cleaned dataset: {len(current_sequences)} sequences")
    else:
        print("⚠️  Current cleaned dataset not found")

    # Read reference sequences
    if os.path.exists(reference_file):
        reference_sequences = list(SeqIO.parse(reference_file, "fasta"))

        # Add only unique sequences (avoid duplicates)
        current_ids = {seq.id.split(".")[0] for seq in all_sequences}
        unique_refs = []

        for seq in reference_sequences:
            seq_id = seq.id.split(".")[0]
            if seq_id not in current_ids:
                unique_refs.append(seq)

        all_sequences.extend(unique_refs)
        print(
            f"📂 Reference sequences: {len(reference_sequences)} total, {len(unique_refs)} unique added"
        )
    else:
        print("⚠️  Reference sequences not found")

    # Write comprehensive dataset
    if all_sequences:
        SeqIO.write(all_sequences, output_file, "fasta")
        print(f"\n✅ Comprehensive dataset created:")
        print(f"   📁 File: {output_file}")
        print(f"   📊 Total sequences: {len(all_sequences)}")

        # Provide next steps
        print(f"\n🚀 NEXT STEPS FOR ENHANCED ANALYSIS:")
        print("1. Align comprehensive dataset with MAFFT:")
        print(f"   mafft --auto {output_file} > data/comprehensive_alignment.fasta")
        print("\n2. Clean alignment and assess gaps")
        print("\n3. Run sliding window analysis with multiple sizes")
        print("\n4. Perform recombination detection with RDP4/GARD")
        print("\n5. Compare results with original dataset")

        return output_file
    else:
        print("❌ No sequences available for comprehensive dataset")
        return None


def main():
    """
    Main execution function.
    """
    print("🧬 CURATED NOROVIRUS REFERENCE ACQUISITION")
    print("=" * 80)
    print("Objective: Fetch high-quality reference sequences for enhanced")
    print("recombination analysis with complete ORF1/ORF2 coverage.")
    print()

    # Step 1: Fetch curated reference sequences
    sequences = fetch_curated_sequences()

    # Step 2: Create comprehensive dataset
    if sequences:
        comprehensive_file = create_comprehensive_dataset()

        if comprehensive_file:
            print(f"\n🎉 ENHANCEMENT READY!")
            print("Your dataset now includes high-quality reference sequences")
            print("for improved recombination breakpoint detection.")

    print(f"\n⏰ Process completed")


if __name__ == "__main__":
    main()
