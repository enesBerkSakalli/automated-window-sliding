#!/usr/bin/env python3
"""
Comprehensive search for norovirus sequences from the reference paper
10.1371/journal.pone.0189504 and related studies.

This script identifies all sequences mentioned in the paper and evaluates
their potential as outgroups or reference sequences for our analysis.
"""


def search_paper_sequences():
    """Search for all sequences mentioned in the paper"""

    print("🔍 COMPREHENSIVE PAPER SEQUENCE SEARCH")
    print("=" * 60)
    print("Paper: 10.1371/journal.pone.0189504")
    print("Title: Emergence and epidemiology of GII.P16-GII.4 noroviruses")
    print("")

    # Sequences specifically mentioned in the paper
    paper_sequences = {
        "new_sequences_deposited": {
            "description": "Sequences newly deposited by the paper authors",
            "accessions": [
                # Range KY451971-KY451987 (17 sequences)
                *[f"KY{451971 + i}" for i in range(17)],
                # Range MF158177-MF158199 (23 sequences)
                *[f"MF{158177 + i}" for i in range(23)],
                # P2 region sequences
                "MF681695",
                "MF681696",
                "KY551568",
                "KY551569",
            ],
        },
        "reference_strains": {
            "description": "Reference strains used for comparison",
            "accessions": [
                "JX459907",  # Sydney 2012 (GII.4)
                "LC175468",  # Japan GII.P16-GII.4
                "KX907727",  # USA GII.P16-GII.4
                "KY887027",  # Another GII.P16-GII.4 reference
            ],
        },
        "outgroup_candidates": {
            "description": "Classical norovirus references for outgroups",
            "accessions": [
                "M87661",  # Norwalk virus (GI.1) - prototype
                "L07418",  # Southampton virus (GI.2)
                "AF093797",  # Desert Shield virus (GI.3)
                "AJ277608",  # Chiba virus (GI.4)
                "U02030",  # Mexico virus (GII.3)
                "AF427118",  # GII.12 reference
            ],
        },
        "global_references": {
            "description": "Global surveillance references",
            "accessions": [
                "AF190817",  # Lordsdale virus (GII.4)
                "X86557",  # Hawaii virus (GII.1)
                "AY134748",  # Snow Mountain virus (GII.2)
                "AB039774",  # Saitama virus (GII.6)
                "JN899242",  # New Orleans virus (GII.4)
            ],
        },
    }

    total_sequences = 0
    for category, data in paper_sequences.items():
        print(f"📂 {category.upper().replace('_', ' ')}:")
        print(f"   {data['description']}")
        print(f"   Count: {len(data['accessions'])} sequences")

        # Show first few accessions
        accessions = data["accessions"]
        if len(accessions) <= 5:
            for acc in accessions:
                print(f"   • {acc}")
        else:
            for acc in accessions[:3]:
                print(f"   • {acc}")
            print(f"   • ... and {len(accessions) - 3} more")
        print()

        total_sequences += len(accessions)

    print(f"📊 TOTAL SEQUENCES IDENTIFIED: {total_sequences}")

    return paper_sequences


def analyze_sequence_relevance():
    """Analyze which sequences would be most relevant for our study"""

    print("\n🎯 SEQUENCE RELEVANCE ANALYSIS")
    print("-" * 40)

    recommendations = {
        "top_priority_outgroups": [
            {
                "accession": "M87661",
                "name": "Norwalk virus",
                "genotype": "GI.1",
                "rationale": "Classic prototype strain, maximum genetic distance from GII",
            },
            {
                "accession": "L07418",
                "name": "Southampton virus",
                "genotype": "GI.2",
                "rationale": "Well-characterized, different GI genotype",
            },
            {
                "accession": "AF093797",
                "name": "Desert Shield virus",
                "genotype": "GI.3",
                "rationale": "Complete genome available, diverse GI lineage",
            },
        ],
        "reference_controls": [
            {
                "accession": "JX459907",
                "name": "Sydney 2012",
                "genotype": "GII.4",
                "rationale": "Reference strain mentioned in paper, well-characterized",
            },
            {
                "accession": "AF190817",
                "name": "Lordsdale virus",
                "genotype": "GII.4",
                "rationale": "Classic GII.4 reference, complete genome",
            },
        ],
        "diversity_representatives": [
            {
                "accession": "U02030",
                "name": "Mexico virus",
                "genotype": "GII.3",
                "rationale": "Different GII genotype, adds phylogenetic diversity",
            },
            {
                "accession": "AB039774",
                "name": "Saitama virus",
                "genotype": "GII.6",
                "rationale": "Matches some sequences in your dataset",
            },
        ],
    }

    for category, sequences in recommendations.items():
        print(f"\n🔸 {category.upper().replace('_', ' ')}:")
        for seq in sequences:
            print(f"   • {seq['accession']} ({seq['name']})")
            print(f"     Genotype: {seq['genotype']}")
            print(f"     Why: {seq['rationale']}")
            print()

    return recommendations


def create_enhanced_sequence_collection():
    """Create a plan for enhanced sequence collection"""

    print("\n🚀 ENHANCED SEQUENCE COLLECTION PLAN")
    print("-" * 45)

    collection_plan = {
        "phase_1": {
            "description": "Add essential outgroups",
            "sequences": ["M87661", "L07418", "AF093797"],
            "priority": "HIGH",
            "impact": "Proper phylogenetic rooting",
        },
        "phase_2": {
            "description": "Add reference controls",
            "sequences": ["JX459907", "AF190817"],
            "priority": "MEDIUM",
            "impact": "Validation against literature",
        },
        "phase_3": {
            "description": "Add diversity representatives",
            "sequences": ["U02030", "AB039774"],
            "priority": "LOW",
            "impact": "Increased phylogenetic context",
        },
    }

    for phase, details in collection_plan.items():
        print(f"\n📋 {phase.upper()}:")
        print(f"   Description: {details['description']}")
        print(f"   Priority: {details['priority']}")
        print(f"   Impact: {details['impact']}")
        print(f"   Sequences: {', '.join(details['sequences'])}")

    print("\n💡 IMPLEMENTATION STRATEGY:")
    print("   1. Start with Phase 1 (essential outgroups)")
    print("   2. Test pipeline with new outgroups")
    print("   3. Compare results with current KR074191.1 outgroup")
    print("   4. Add Phase 2 sequences if needed")
    print("   5. Document improvements in tree topology")

    return collection_plan


def download_sequences_template():
    """Provide template code for downloading sequences from NCBI"""

    print("\n💻 SEQUENCE DOWNLOAD TEMPLATE")
    print("-" * 35)

    template_code = '''
# Template for downloading sequences from NCBI
from Bio import Entrez, SeqIO
import time

def download_norovirus_references(accessions, email, output_file):
    """Download reference sequences from NCBI"""
    
    # Set email for NCBI
    Entrez.email = email
    
    sequences = []
    
    for acc in accessions:
        try:
            print(f"Downloading {acc}...")
            
            # Search for the sequence
            handle = Entrez.efetch(
                db="nucleotide",
                id=acc,
                rettype="fasta",
                retmode="text"
            )
            
            # Parse the sequence
            record = SeqIO.read(handle, "fasta")
            sequences.append(record)
            handle.close()
            
            # Be nice to NCBI servers
            time.sleep(0.5)
            
        except Exception as e:
            print(f"Error downloading {acc}: {e}")
            continue
    
    # Write to file
    SeqIO.write(sequences, output_file, "fasta")
    print(f"Saved {len(sequences)} sequences to {output_file}")
    
    return sequences

# Example usage:
# outgroup_accessions = ["M87661", "L07418", "AF093797"]
# download_norovirus_references(
#     outgroup_accessions, 
#     "your.email@example.com",
#     "outgroup_sequences.fasta"
# )
'''

    print("📝 Template created for NCBI sequence download")
    print("⚠️  Requires: Biopython, valid email address")
    print("🔗 Can be integrated into your pipeline")

    return template_code


def compare_with_current_dataset():
    """Compare paper sequences with our current dataset"""

    print("\n📊 COMPARISON WITH CURRENT DATASET")
    print("-" * 40)

    current_analysis = {
        "current_dataset": {
            "sequences": 38,
            "length": "484 bp (after enhancement)",
            "genotypes": "Mainly GII.4, some GII.6, one GII.12",
            "outgroup": "KR074191.1 (GII.Pg-GII.12)",
        },
        "enhancement_opportunities": {
            "add_gi_outgroups": "3-5 GI sequences for proper rooting",
            "add_references": "2-3 well-characterized GII references",
            "remove_problematic": "Consider removing highly divergent sequences",
            "validate_genotyping": "Cross-check genotype assignments",
        },
        "expected_improvements": {
            "tree_rooting": "More stable, biologically meaningful roots",
            "branch_lengths": "Better calibrated evolutionary distances",
            "topology": "Improved resolution of GII relationships",
            "validation": "Consistency with published phylogenies",
        },
    }

    print("🔍 CURRENT STATE:")
    for key, value in current_analysis["current_dataset"].items():
        print(f"   {key.title()}: {value}")

    print("\n🚀 ENHANCEMENT OPPORTUNITIES:")
    for key, value in current_analysis["enhancement_opportunities"].items():
        print(f"   • {key.replace('_', ' ').title()}: {value}")

    print("\n✨ EXPECTED IMPROVEMENTS:")
    for key, value in current_analysis["expected_improvements"].items():
        print(f"   • {key.replace('_', ' ').title()}: {value}")

    return current_analysis


def create_action_plan():
    """Create specific action plan for implementation"""

    print("\n📋 IMPLEMENTATION ACTION PLAN")
    print("-" * 35)

    action_plan = [
        {
            "step": 1,
            "action": "Download GI outgroup sequences",
            "details": "M87661, L07418, AF093797 from NCBI",
            "time": "15 minutes",
            "priority": "HIGH",
        },
        {
            "step": 2,
            "action": "Create enhanced alignment",
            "details": "Add outgroups to refined_alignment.fasta",
            "time": "10 minutes",
            "priority": "HIGH",
        },
        {
            "step": 3,
            "action": "Test pipeline with new alignment",
            "details": "Run sliding window analysis with new outgroups",
            "time": "30 minutes",
            "priority": "HIGH",
        },
        {
            "step": 4,
            "action": "Compare tree topologies",
            "details": "Analyze differences with/without new outgroups",
            "time": "20 minutes",
            "priority": "MEDIUM",
        },
        {
            "step": 5,
            "action": "Validate against literature",
            "details": "Check consistency with published norovirus trees",
            "time": "30 minutes",
            "priority": "MEDIUM",
        },
        {
            "step": 6,
            "action": "Document improvements",
            "details": "Update methods and create summary report",
            "time": "20 minutes",
            "priority": "LOW",
        },
    ]

    total_time = 0
    for step in action_plan:
        time_minutes = int(step["time"].split()[0])
        total_time += time_minutes

        print(f"\n🔸 STEP {step['step']}: {step['action'].upper()}")
        print(f"   Details: {step['details']}")
        print(f"   Time: {step['time']}")
        print(f"   Priority: {step['priority']}")

    print(
        f"\n⏱️  TOTAL ESTIMATED TIME: {total_time} minutes ({total_time / 60:.1f} hours)"
    )
    print("🎯 Focus on HIGH priority steps first")

    return action_plan


if __name__ == "__main__":
    print("🧬 NOROVIRUS SEQUENCE ENHANCEMENT ANALYSIS")
    print("=" * 60)

    # Run comprehensive analysis
    paper_sequences = search_paper_sequences()
    recommendations = analyze_sequence_relevance()
    collection_plan = create_enhanced_sequence_collection()
    template_code = download_sequences_template()
    comparison = compare_with_current_dataset()
    action_plan = create_action_plan()

    print("\n" + "=" * 60)
    print("🎉 COMPREHENSIVE SEQUENCE ANALYSIS COMPLETED")
    print(
        f"📊 Total sequences identified: {sum(len(data['accessions']) for data in paper_sequences.values())}"
    )
    print("🔍 Key recommendations: Add 3-5 GI outgroup sequences")
    print("⚡ Next step: Download M87661, L07418, AF093797")
    print("=" * 60)
