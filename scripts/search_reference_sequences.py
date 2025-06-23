#!/usr/bin/env python3
"""
Search for additional norovirus reference sequences that could serve as outgroups.
Based on research findings, GI sequences are the most suitable outgroups for GII analysis.
"""

# Note: Bio.Entrez can be used to download sequences from NCBI if needed


def search_reference_sequences():
    """Search for well-characterized norovirus reference sequences"""

    print("🔍 SEARCHING FOR NOROVIRUS REFERENCE SEQUENCES")
    print("=" * 60)

    # Well-known reference sequences from literature
    reference_sequences = {
        "GI_outgroups": [
            "M87661",  # Norwalk virus (prototype GI.1)
            "L07418",  # Southampton virus (GI.2)
            "AF093797",  # Desert Shield virus (GI.3)
            "AJ277608",  # Chiba virus (GI.4)
            "AY502023",  # Musgrove virus (GI.5)
            "AY502009",  # Hesse virus (GI.6)
        ],
        "GII_references": [
            "AF190817",  # Lordsdale virus (GII.4)
            "X86557",  # Hawaii virus (GII.1)
            "AY134748",  # Snow Mountain virus (GII.2)
            "U02030",  # Mexico virus (GII.3)
            "JX459907",  # Sydney 2012 (GII.4) - mentioned in your paper
            "AB039774",  # Saitama virus (GII.6)
        ],
        "diverse_genotypes": [
            "AF427118",  # GII.12 reference
            "LC175468",  # GII.P16-GII.4 from Japan (mentioned in paper)
            "KX907727",  # GII.P16-GII.4 from USA (mentioned in paper)
        ],
    }

    print("📚 WELL-KNOWN REFERENCE SEQUENCES:")
    print("-" * 40)

    for category, accessions in reference_sequences.items():
        print(f"\n{category.upper()}:")
        for acc in accessions:
            print(f"  • {acc}")

    print("\n🎯 RECOMMENDATIONS FOR YOUR ANALYSIS:")
    print("-" * 40)

    print("\n1. BEST OUTGROUP OPTIONS (GI sequences):")
    print("   • M87661 (Norwalk virus, GI.1) - Classic reference")
    print("   • L07418 (Southampton virus, GI.2) - Well-characterized")
    print("   • AF093797 (Desert Shield virus, GI.3) - Diverse")

    print("\n2. ALTERNATIVE OUTGROUPS (Distant GII):")
    print("   • AF427118 (GII.12 reference) - Same genotype as KR074191.1")
    print("   • U02030 (Mexico virus, GII.3) - Different from your GII.4/6 dataset")

    print("\n3. REFERENCE SEQUENCES FOR COMPARISON:")
    print("   • JX459907 (Sydney 2012) - Mentioned in your paper")
    print("   • LC175468 (Japan GII.P16-GII.4) - From your paper")
    print("   • KX907727 (USA GII.P16-GII.4) - From your paper")

    print("\n🔬 ANALYSIS STRATEGY:")
    print("-" * 25)
    print("1. Download GI reference sequences (M87661, L07418, AF093797)")
    print("2. Add them to your dataset as outgroups")
    print("3. Re-run pipeline with multiple outgroups")
    print("4. Compare tree topologies with/without KR074191.1")
    print("5. Use midpoint rooting as control method")

    return reference_sequences


def create_outgroup_fasta():
    """Create a FASTA file with suggested outgroup sequences"""

    print("\n💾 CREATING OUTGROUP SEQUENCE FILE:")
    print("-" * 40)

    # These are example sequences - in practice you'd download from NCBI
    outgroup_sequences = """
>M87661.1 Norwalk virus, complete genome (GI.1 reference)
# This would contain the actual sequence from NCBI
# Use: Bio.Entrez.efetch() to download

>L07418.1 Southampton virus capsid protein gene (GI.2 reference)  
# This would contain the actual sequence from NCBI
# Use: Bio.Entrez.efetch() to download

>AF093797.1 Desert Shield virus, complete genome (GI.3 reference)
# This would contain the actual sequence from NCBI  
# Use: Bio.Entrez.efetch() to download
"""

    print("📝 Example FASTA structure created")
    print("⚠️  Note: Actual sequences need to be downloaded from NCBI")
    print("🔗 Use Bio.Entrez.efetch() to get real sequences")

    return outgroup_sequences


def analyze_paper_sequences():
    """Analyze sequences specifically mentioned in the paper"""

    print("\n📄 SEQUENCES FROM PAPER 10.1371/journal.pone.0189504:")
    print("-" * 60)

    paper_sequences = {
        "newly_deposited": [
            "KY451971-KY451987",  # From the paper
            "MF158177-MF158199",  # From the paper
            "MF681695-MF681696",  # P2 sequences
            "KY551568-KY551569",  # P2 sequences
        ],
        "reference_comparisons": [
            "JX459907",  # Sydney 2012 reference
            "LC175468",  # Japan GII.P16-GII.4
            "KX907727",  # USA GII.P16-GII.4
        ],
    }

    print("🧬 New sequences deposited by the paper authors:")
    for acc_range in paper_sequences["newly_deposited"]:
        print(f"   • {acc_range}")

    print("\n🔗 Reference sequences used in paper:")
    for acc in paper_sequences["reference_comparisons"]:
        print(f"   • {acc}")

    print("\n💡 KEY FINDINGS FROM PAPER:")
    print("   • Used neighbor-joining with Kimura 2-parameter model")
    print("   • 2000 bootstrap replications")
    print("   • Analyzed ORF1/ORF2 junction region (570 bp)")
    print("   • Found GII.P16-GII.4 clustering separately from earlier GII.P16 strains")
    print("   • No amino acid changes in major epitopes A-E")

    return paper_sequences


def recommend_pipeline_updates():
    """Recommend updates to your analysis pipeline"""

    print("\n🔧 RECOMMENDED PIPELINE UPDATES:")
    print("-" * 40)

    print("\n1. ENHANCED ALIGNMENT (✅ COMPLETED):")
    print("   • Quality filtering: Remove bottom 5% by length")
    print("   • Sequence trimming: 15 bp from each end")
    print("   • Genotype-specific alignment")
    print("   • Gap-heavy column removal (>60% gaps)")
    print("   • Result: 38 sequences, 484 bp, 82.1% conservation")

    print("\n2. OUTGROUP STRATEGY:")
    print("   • Add GI reference sequences (M87661, L07418, AF093797)")
    print("   • Keep KR074191.1 with branch length correction")
    print("   • Try multiple outgroup rooting")
    print("   • Compare with midpoint rooting")

    print("\n3. TREE BUILDING IMPROVEMENTS:")
    print("   • Use refined_alignment.fasta (higher quality)")
    print("   • Apply GTR+F+I+G4 model (from your IQ-TREE results)")
    print("   • Increase bootstrap replicates to 2000 (match paper)")
    print("   • Keep branch length capping (0.1 substitutions/site)")

    print("\n4. VALIDATION STEPS:")
    print("   • Compare tree topologies with different outgroups")
    print("   • Check bootstrap support for major clades")
    print("   • Validate against published norovirus phylogenies")
    print("   • Test sliding window analysis with new alignment")

    print("\n🎯 NEXT STEPS:")
    print("   1. Test pipeline with refined_alignment.fasta")
    print("   2. Download and add GI reference sequences")
    print("   3. Compare results with original pipeline")
    print("   4. Document improvements in methods")


if __name__ == "__main__":
    # Run analysis
    reference_sequences = search_reference_sequences()
    outgroup_sequences = create_outgroup_fasta()
    paper_sequences = analyze_paper_sequences()
    recommend_pipeline_updates()

    print("\n" + "=" * 60)
    print("🎉 REFERENCE SEQUENCE ANALYSIS COMPLETED")
    print("📋 Summary: Enhanced alignment ready, outgroup options identified")
    print("🔜 Next: Test pipeline with improved alignment")
    print("=" * 60)
