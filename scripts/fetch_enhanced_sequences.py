#!/usr/bin/env python3
"""
Enhanced Norovirus Sequence Acquisition Script
Fetches high-quality, recent complete/near-complete genome sequences
for improved recombination analysis.

Based on research findings:
- Focus on complete genomes with full ORF1/ORF2 coverage
- Include recent recombinant strains (GII.P16, GII.P21, GII.P31)
- Target sequences from 2020-2024 for current relevance
- Prioritize sequences with known recombination breakpoints
"""

from Bio import Entrez, SeqIO
import time
import os
from datetime import datetime

# Set your email for NCBI Entrez
Entrez.email = "your.email@example.com"


def fetch_recent_complete_genomes():
    """
    Fetch recent complete/near-complete norovirus genomes from GenBank.
    Focus on high-quality sequences for recombination analysis.
    """

    # High-priority recent complete genome accessions identified from research
    priority_accessions = [
        # Recent complete genomes (2020-2024)
        "MK753033",  # Referenced in Illumina-Nanopore study
        "KY210980",  # Complete GII.P16-GII.4/RUS/Novosibirsk/NS16-C38/2016
        "KY421121",  # Complete GII.P16/GII.2 from China 2016
        "MK789430",  # Complete genome from Valencia, Spain
        "MK789431",  # Complete genome from Valencia, Spain
        "MK789432",  # Complete genome from Valencia, Spain
        "MK789433",  # Complete genome from Valencia, Spain
        # Additional high-quality references
        "KY210977",  # GII.Pe-GII.4/RUS/Novosibirsk/NS16-C13/2016
        "KY210919",  # GII.P21-GII.3/RUS/Novosibirsk/NS16-C23/2016
    ]

    # Search terms for additional recent sequences
    search_terms = [
        # Recent complete genomes
        '("complete genome"[Title] OR "near complete genome"[Title]) AND norovirus[Title] AND 2020:2024[Publication Date]',
        # Recombinant strains
        'norovirus[Title] AND (GII.P16[Title] OR GII.P21[Title] OR GII.P31[Title]) AND "complete"[Title] AND 2018:2024[Publication Date]',
        # High-quality ORF1/ORF2 sequences
        "norovirus[Title] AND (ORF1[Title] AND ORF2[Title]) AND 2020:2024[Publication Date]",
        # Recent recombination studies
        'norovirus[Title] AND recombinant[Title] AND "complete genome"[Title] AND 2019:2024[Publication Date]',
    ]

    print("🔍 SEARCHING FOR ENHANCED NOROVIRUS SEQUENCES")
    print("=" * 60)

    all_accessions = set(priority_accessions)

    # Search for additional sequences
    for i, search_term in enumerate(search_terms, 1):
        print(f"\n🔎 Search {i}: {search_term[:80]}...")
        try:
            handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=50)
            search_results = Entrez.read(handle)
            handle.close()

            new_ids = search_results["IdList"]
            print(f"   Found {len(new_ids)} sequences")

            # Get accession numbers
            if new_ids:
                handle = Entrez.efetch(db="nucleotide", id=new_ids, rettype="acc")
                accessions = handle.read().strip().split("\n")
                handle.close()

                all_accessions.update(accessions)
                time.sleep(1)  # Be nice to NCBI

        except Exception as e:
            print(f"   Error in search {i}: {e}")
            continue

    print(f"\n📊 Total unique accessions found: {len(all_accessions)}")
    return list(all_accessions)


def fetch_sequences_batch(accessions, output_file, batch_size=20):
    """
    Fetch sequences in batches to avoid overwhelming NCBI servers.
    """
    print(f"\n💾 FETCHING SEQUENCES TO {output_file}")
    print("=" * 60)

    sequences = []
    failed_accessions = []

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        print(f"\n📦 Processing batch {i // batch_size + 1}: {len(batch)} sequences")

        try:
            # Fetch sequences
            handle = Entrez.efetch(
                db="nucleotide", id=batch, rettype="fasta", retmode="text"
            )

            batch_sequences = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            # Filter for high-quality sequences
            for seq in batch_sequences:
                seq_len = len(seq.seq)

                # Filter criteria for norovirus complete/near-complete genomes
                if seq_len >= 6000:  # Minimum length for near-complete genomes
                    sequences.append(seq)
                    print(f"   ✅ {seq.id}: {seq_len}bp - {seq.description[:60]}...")
                else:
                    print(f"   ❌ {seq.id}: {seq_len}bp - Too short")
                    failed_accessions.append(seq.id)

            time.sleep(2)  # Be respectful to NCBI

        except Exception as e:
            print(f"   ❌ Batch failed: {e}")
            failed_accessions.extend(batch)

    # Write sequences to file
    if sequences:
        SeqIO.write(sequences, output_file, "fasta")
        print(
            f"\n✅ SUCCESS: {len(sequences)} high-quality sequences saved to {output_file}"
        )
    else:
        print("\n❌ No sequences retrieved")

    if failed_accessions:
        print(f"\n⚠️  Failed to retrieve {len(failed_accessions)} accessions:")
        for acc in failed_accessions[:10]:  # Show first 10
            print(f"   - {acc}")
        if len(failed_accessions) > 10:
            print(f"   ... and {len(failed_accessions) - 10} more")

    return sequences, failed_accessions


def analyze_sequence_quality(sequences):
    """
    Analyze the quality and characteristics of retrieved sequences.
    """
    print(f"\n📈 SEQUENCE QUALITY ANALYSIS")
    print("=" * 60)

    if not sequences:
        print("No sequences to analyze")
        return

    lengths = [len(seq.seq) for seq in sequences]
    genogroups = {}

    for seq in sequences:
        header = seq.description.upper()

        # Extract genogroup information
        if "GII.P" in header:
            import re

            geno_match = re.search(r"GII\.P\d+[A-Z]*", header)
            if geno_match:
                geno = geno_match.group()
                genogroups[geno] = genogroups.get(geno, 0) + 1

    print(f"📊 Sequence Statistics:")
    print(f"   Total sequences: {len(sequences)}")
    print(f"   Length range: {min(lengths):,} - {max(lengths):,} bp")
    print(f"   Average length: {sum(lengths) // len(lengths):,} bp")
    print(f"   Complete genomes (>7000bp): {sum(1 for l in lengths if l > 7000)}")
    print(
        f"   Near-complete (6000-7000bp): {sum(1 for l in lengths if 6000 <= l <= 7000)}"
    )

    if genogroups:
        print(f"\n🧬 Genogroup Distribution:")
        for geno, count in sorted(genogroups.items()):
            print(f"   {geno}: {count} sequences")


def create_enhanced_dataset():
    """
    Create an enhanced dataset by combining current sequences with new high-quality sequences.
    """
    print(f"\n🔧 CREATING ENHANCED DATASET")
    print("=" * 60)

    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    current_file = os.path.join(base_dir, "data", "cleaned_alignment_combined.fasta")
    new_sequences_file = os.path.join(
        base_dir, "data", "enhanced_complete_genomes.fasta"
    )
    final_combined_file = os.path.join(base_dir, "data", "ultra_enhanced_dataset.fasta")

    # Read current cleaned alignment
    if os.path.exists(current_file):
        current_sequences = list(SeqIO.parse(current_file, "fasta"))
        print(f"📂 Current dataset: {len(current_sequences)} sequences")
    else:
        print("❌ Current cleaned alignment not found")
        return

    # Read new sequences
    if os.path.exists(new_sequences_file):
        new_sequences = list(SeqIO.parse(new_sequences_file, "fasta"))
        print(f"📂 New sequences: {len(new_sequences)} sequences")
    else:
        print("❌ New sequences file not found")
        return

    # Combine datasets (avoid duplicates by accession)
    current_ids = {seq.id.split(".")[0] for seq in current_sequences}
    unique_new_sequences = []

    for seq in new_sequences:
        seq_id = seq.id.split(".")[0]
        if seq_id not in current_ids:
            unique_new_sequences.append(seq)

    print(f"🔄 Adding {len(unique_new_sequences)} unique new sequences")

    # Combine all sequences
    all_sequences = current_sequences + unique_new_sequences

    # Write combined dataset
    SeqIO.write(all_sequences, final_combined_file, "fasta")

    print(f"✅ Enhanced dataset created: {final_combined_file}")
    print(f"   Total sequences: {len(all_sequences)}")
    print(f"   Original: {len(current_sequences)}")
    print(f"   Added: {len(unique_new_sequences)}")

    return final_combined_file


def main():
    """
    Main function to execute enhanced sequence acquisition.
    """
    print("🧬 ENHANCED NOROVIRUS SEQUENCE ACQUISITION")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("\nObjective: Acquire high-quality complete/near-complete norovirus genomes")
    print("for enhanced recombination analysis with better ORF1/ORF2 coverage.")

    # Create output directory
    output_dir = "/Users/berksakalli/Projects/automated-window-sliding/data"
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Search for sequences
    accessions = fetch_recent_complete_genomes()

    if not accessions:
        print("❌ No accessions found")
        return

    # Step 2: Fetch sequences
    output_file = os.path.join(output_dir, "enhanced_complete_genomes.fasta")
    sequences, failed = fetch_sequences_batch(accessions, output_file)

    # Step 3: Analyze quality
    analyze_sequence_quality(sequences)

    # Step 4: Create enhanced dataset
    if sequences:
        enhanced_file = create_enhanced_dataset()

        print(f"\n🎉 ENHANCEMENT COMPLETE!")
        print("=" * 60)
        print("Next steps:")
        print("1. Align the enhanced dataset with MAFFT")
        print("2. Clean alignment and remove problematic sequences")
        print("3. Run sliding window analysis with multiple window sizes")
        print("4. Perform recombination analysis with RDP4/GARD")
        print(f"\nEnhanced dataset ready: {enhanced_file}")

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
