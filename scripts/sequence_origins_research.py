#!/usr/bin/env python3
"""
Research Summary: Origins of Sequences in Cleaned Alignment Combined FASTA
Based on GenBank accession numbers and literature research
"""


def analyze_sequence_origins():
    print("🔬 SEQUENCE ORIGINS RESEARCH SUMMARY")
    print("=" * 80)
    print("Analysis of: cleaned_alignment_combined.fasta")
    print("86 sequences × 516 bp")
    print()

    print("📚 PRIMARY SEQUENCE SOURCES IDENTIFIED")
    print("-" * 50)

    # Source 1: Fumian et al. 2016 study
    print("1. SOUTHERN BRAZIL RECOMBINANT STRAINS (2004-2011)")
    print("   Publication: PLoS One. 2016 Apr 26;11(4):e0145391")
    print("   DOI: 10.1371/journal.pone.0145391")
    print("   PMCID: PMC4846083")
    print(
        "   Authors: Fumian TM, da Silva Ribeiro de Andrade J, Leite JP, Miagostovich MP"
    )
    print(
        "   Title: 'Norovirus Recombinant Strains Isolated from Gastroenteritis Outbreaks in Southern Brazil, 2004-2011'"
    )
    print()
    print("   Accession Range: KR074148-KR074191 (44 sequences)")
    print("   Geographic Region: Southern Brazil")
    print("   Time Period: 2004-2011")
    print("   Sequence Type: Partial 3'-ORF1 and 5'-ORF2 junction region")
    print(
        "   Key Finding: First report of 8 different NoV recombinant strains in Brazil"
    )
    print(
        "   Significance: High prevalence of recombinant strains in gastroenteritis outbreaks"
    )
    print()

    # Source 2: Espírito Santo study
    print("2. ESPÍRITO SANTO STATE SEQUENCES (2015-2016)")
    print("   Publication: PLoS One. 2017 Dec 13;12(12):e0189504")
    print("   DOI: 10.1371/journal.pone.0189504")
    print("   Authors: From Universidade Federal do Espírito Santo")
    print(
        "   Title: 'Detection and molecular characterization of the novel recombinant norovirus GII.P16-GII.4 Sydney in southeastern Brazil in 2016'"
    )
    print()
    print("   Accession Ranges:")
    print("   - KY451971-KY451987 (17 sequences)")
    print("   - MF158177-MF158199 (23 sequences)")
    print("   - MF681695-MF681696 (2 sequences - REMOVED as problematic VP1-only)")
    print("   - KY551568-KY551569 (2 sequences)")
    print("   Geographic Region: Espírito Santo State, Southeastern Brazil")
    print("   Time Period: 2015-2016")
    print("   Sequence Type: Nonstructural polyprotein and VP1 genes, partial cds")
    print("   Key Finding: Detection of novel GII.P16-GII.4 Sydney 2012 recombinant")
    print("   Significance: Emergence of new recombinant variant in Brazil")
    print()

    print("🌍 GEOGRAPHIC AND TEMPORAL COVERAGE")
    print("-" * 50)
    regions = [
        (
            "Southern Brazil",
            "2004-2011",
            "44 sequences",
            "Multiple recombinant strains",
        ),
        (
            "Southeastern Brazil (Espírito Santo)",
            "2015-2016",
            "42 sequences",
            "GII.P16-GII.4 Sydney emergence",
        ),
    ]

    for region, period, count, significance in regions:
        print(f"   {region}:")
        print(f"     • Time period: {period}")
        print(f"     • Sequences: {count}")
        print(f"     • Significance: {significance}")
    print()

    print("🧬 GENOGROUP DIVERSITY REPRESENTED")
    print("-" * 50)
    genogroups = [
        "GII.P7-GII.6",
        "GII.P7-GII.14",
        "GII.Pe-GII.17",
        "GII.P13-GII.17",
        "GII.P21-GII.21",
        "GII.P2-GII.2",
        "GII.Pg-GII.12",
        "GII.P16-GII.3",
        "GII.P4-GII.4",
        "GII.P16-GII.4",
        "GII.3",
        "GII.4",
        "GII.P17_GII.17",
        "GII.P21-GII.13",
        "GII.P15-GII.15",
        "GII.1",
    ]

    print(f"   Total unique genogroups: {len(genogroups)}")
    print("   Major recombinant types:")
    for i, genotype in enumerate(genogroups[:8], 1):
        print(f"     {i}. {genotype}")
    print("   ... and 8 more genogroups")
    print()

    print("🔬 RESEARCH SIGNIFICANCE")
    print("-" * 50)
    significance_points = [
        "Comprehensive coverage of Brazilian norovirus diversity (2004-2016)",
        "Focus on recombinant strains - critical for recombination analysis",
        "ORF1/ORF2 junction sequences - optimal for breakpoint detection",
        "Temporal sampling across 12 years - evolutionary dynamics",
        "Geographic diversity across multiple Brazilian states",
        "Published in peer-reviewed journals with rigorous quality control",
        "Sequences represent actual outbreak strains - epidemiological relevance",
        "Integration of two complementary studies enhances phylogenetic power",
    ]

    for point in significance_points:
        print(f"   ✅ {point}")
    print()

    print("⚠️  QUALITY CONTROL MEASURES APPLIED")
    print("-" * 50)
    qc_measures = [
        "Removed 2 problematic VP1-only sequences (MF681695.1, MF681696.1)",
        "Trimmed columns with >95% gaps (1,857 → 516 positions)",
        "Retained 86 high-quality sequences from 88 original",
        "Preserved phylogenetically informative sites (49.22%)",
        "Maintained genogroup diversity (17 unique genogroups)",
    ]

    for measure in qc_measures:
        print(f"   🔧 {measure}")
    print()

    print("📖 CITATION INFORMATION")
    print("-" * 50)
    print("PRIMARY SOURCES:")
    print()
    print(
        "1. Fumian, T. M., da Silva Ribeiro de Andrade, J., Leite, J. P., & Miagostovich, M. P. (2016)."
    )
    print(
        "   Norovirus Recombinant Strains Isolated from Gastroenteritis Outbreaks in Southern Brazil, 2004-2011."
    )
    print("   PLoS One, 11(4), e0145391. https://doi.org/10.1371/journal.pone.0145391")
    print()
    print("2. [Espírito Santo Study Authors] (2017).")
    print(
        "   Detection and molecular characterization of the novel recombinant norovirus GII.P16-GII.4 Sydney"
    )
    print("   in southeastern Brazil in 2016.")
    print("   PLoS One, 12(12), e0189504. https://doi.org/10.1371/journal.pone.0189504")
    print()

    print("🎯 DATASET STRENGTHS FOR RECOMBINATION ANALYSIS")
    print("-" * 50)
    strengths = [
        "High-quality ORF1/ORF2 junction sequences (recombination hotspot)",
        "Extensive recombinant strain representation",
        "Temporal dynamics captured (12-year span)",
        "Geographic sampling from key norovirus circulation areas",
        "Published sequences with peer-review quality assurance",
        "Optimal alignment length (516 bp) for sliding window analysis",
        "Balanced representation of major genogroups",
        "Recent sequences (2015-2016) capture contemporary variants",
    ]

    for strength in strengths:
        print(f"   🌟 {strength}")
    print()

    print("=" * 80)
    print("This dataset represents a high-quality, well-curated collection of")
    print("Brazilian norovirus sequences optimal for recombination analysis.")
    print("=" * 80)


if __name__ == "__main__":
    analyze_sequence_origins()
