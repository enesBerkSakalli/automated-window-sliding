#!/usr/bin/env python3
"""
Focused script to apply improved genotype mappings only to specific tree files:
- all_rooted_trees.newick/.nexus
- rooted_trees_collection.newick/.nexus

This script:
1. Keeps accession numbers (important for identification)
2. Removes duplicates to improve readability
3. Creates better file names for outputs
4. Only processes the specified tree files
"""

import os
import re


def extract_genotype_info(header):
    """Enhanced genotype extraction with better parsing"""
    # Remove the '>' character
    header = header.strip().lstrip(">")

    # Initialize variables
    accession = ""
    genotype = "Unknown"
    p_type = "Unknown"
    year = "Unknown"
    country = "Unknown"
    strain = "Unknown"

    # Extract accession number (first part before space)
    parts = header.split()
    if parts:
        accession = parts[0]

    # Enhanced regex patterns for genotype extraction
    patterns = [
        # Standard format: GII.P4-GII.4, GII.P7-GII.6, etc.
        r"GII\.P(\d+)[_-]GII\.(\d+)",
        # Pe format: GII.Pe-GII.4, GII.Pe_GII.4
        r"GII\.Pe[_-]GII\.(\d+)",
        # Pg format: GII.Pg-GII.12
        r"GII\.Pg[_-]GII\.(\d+)",
        # P13 format: GII.P13-GII.17
        r"GII\.P(\d+)[_-]GII\.(\d+)",
        # Simple GII.4 format
        r"GII\.(\d+)(?:[^P\d]|$)",
        # P16_GII.3 format (underscores)
        r"GII\.P(\d+)_GII\.(\d+)",
        # P17_GII.17 format
        r"P(\d+)_GII\.(\d+)",
    ]

    for pattern in patterns:
        match = re.search(pattern, header)
        if match:
            groups = match.groups()
            if len(groups) == 2:
                p_num, g_num = groups
                p_type = f"GII.P{p_num}"
                genotype = f"GII.{g_num}"
                break
            elif len(groups) == 1:
                g_num = groups[0]
                genotype = f"GII.{g_num}"
                if "Pe" in header:
                    p_type = "GII.Pe"
                elif "Pg" in header:
                    p_type = "GII.Pg"
                break

    # Extract year
    year_match = re.search(r"(\d{4})", header)
    if year_match:
        year = year_match.group(1)

    # Extract country (look for 3-letter country codes)
    country_match = re.search(r"/([A-Z]{2,3})/", header)
    if country_match:
        country = country_match.group(1)

    # Extract strain (last part after /)
    strain_match = re.search(r"/([^/\s]+)\s", header)
    if strain_match:
        strain = strain_match.group(1)

    return {
        "accession": accession,
        "genotype": genotype,
        "p_type": p_type,
        "year": year,
        "country": country,
        "strain": strain,
    }


def create_improved_mapping(fasta_file):
    """Create improved mapping with better display names"""
    accession_to_info = {}
    seen_combinations = set()

    print("🔍 Analyzing FASTA sequences...")

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                info = extract_genotype_info(line)
                accession = info["accession"]

                # Create a comprehensive but readable display name
                if info["genotype"] != "Unknown" and info["p_type"] != "Unknown":
                    # Full dual typing
                    genotype_str = f"{info['p_type']}-{info['genotype']}"
                elif info["genotype"] != "Unknown":
                    # Just genotype
                    genotype_str = info["genotype"]
                elif info["p_type"] != "Unknown":
                    # Just P-type
                    genotype_str = info["p_type"]
                else:
                    genotype_str = "Unknown"

                # Create display name: Accession_Genotype_Year_Country
                display_parts = [accession, genotype_str]
                if info["year"] != "Unknown":
                    display_parts.append(info["year"])
                if info["country"] != "Unknown":
                    display_parts.append(info["country"])

                display_name = "_".join(display_parts)

                # Check for duplicates and make unique if needed
                original_display = display_name
                counter = 1
                while display_name in seen_combinations:
                    display_name = f"{original_display}_{counter}"
                    counter += 1

                seen_combinations.add(display_name)
                accession_to_info[accession] = {
                    "display_name": display_name,
                    "genotype_info": genotype_str,
                    **info,
                }

    print(f"✅ Created mappings for {len(accession_to_info)} sequences")
    return accession_to_info


def process_tree_file(tree_file, accession_mapping, output_suffix=""):
    """Process a tree file with improved mapping"""
    if not os.path.exists(tree_file):
        print(f"⚠️  Tree file not found: {tree_file}")
        return None

    print(f"🌳 Processing: {os.path.basename(tree_file)}")

    with open(tree_file, "r") as f:
        tree_content = f.read()

    # Count replacements for reporting
    replacements_made = 0

    # Replace accession numbers with improved display names
    for accession, info in accession_mapping.items():
        # Look for the accession in various formats
        patterns = [
            accession,  # Just the accession
            f"{accession}:",  # With colon
            f"{accession}_",  # With underscore (if already modified)
        ]

        for pattern in patterns:
            if pattern in tree_content:
                tree_content = tree_content.replace(pattern, info["display_name"])
                if pattern == accession:  # Only count primary replacements
                    replacements_made += 1

    # Create better output filename
    base_dir = os.path.dirname(tree_file)
    base_name = os.path.splitext(os.path.basename(tree_file))[0]
    extension = os.path.splitext(tree_file)[1]

    # Better naming convention
    if "rooted_trees_collection" in base_name:
        output_name = f"rooted_trees_with_genotypes{output_suffix}{extension}"
    elif "all_rooted_trees" in base_name:
        output_name = f"all_rooted_trees_with_genotypes{output_suffix}{extension}"
    else:
        output_name = f"{base_name}_with_genotypes{output_suffix}{extension}"

    output_path = os.path.join(base_dir, output_name)

    # Write the modified tree
    with open(output_path, "w") as f:
        f.write(tree_content)

    print(f"   ✅ Created: {output_name} ({replacements_made} replacements)")
    return output_path


def process_results_directory(results_dir, fasta_file):
    """Process specific tree files in a results directory"""
    print(f"\n🔬 Processing directory: {os.path.basename(results_dir)}")
    print("=" * 60)

    if not os.path.exists(results_dir):
        print(f"❌ Directory not found: {results_dir}")
        return

    # Create improved mapping
    accession_mapping = create_improved_mapping(fasta_file)

    # Target specific files
    target_files = [
        "all_rooted_trees.newick",
        "all_rooted_trees.nexus",
        "rooted_trees_collection.newick",
        "rooted_trees_collection.nexus",
    ]

    processed_files = []

    for target_file in target_files:
        tree_file = os.path.join(results_dir, target_file)
        output_file = process_tree_file(tree_file, accession_mapping)
        if output_file:
            processed_files.append(output_file)

    # Create mapping summary file
    mapping_file = os.path.join(results_dir, "genotype_mapping_summary.txt")
    with open(mapping_file, "w") as f:
        f.write("# Genotype Mapping Summary\n")
        f.write(f"# Generated on: {os.popen('date').read().strip()}\n")
        f.write(f"# Total sequences mapped: {len(accession_mapping)}\n\n")
        f.write("Accession\tDisplay_Name\tGenotype\tP_Type\tYear\tCountry\n")

        for accession, info in sorted(accession_mapping.items()):
            f.write(
                f"{accession}\t{info['display_name']}\t{info['genotype']}\t{info['p_type']}\t{info['year']}\t{info['country']}\n"
            )

    print(f"✅ Created mapping summary: {os.path.basename(mapping_file)}")
    print(
        f"📊 Summary: {len(processed_files)} tree files processed, {len(accession_mapping)} sequences mapped"
    )

    return processed_files


def main():
    """Main function to process all three target directories"""
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    fasta_file = os.path.join(base_dir, "data", "cleaned_alignment_combined.fasta")

    # Target directories
    target_dirs = [
        "results/cleaned_alignment_combined_w250_s10_GTR_F_I_G4_20250619_103929",
        "results/cleaned_alignment_combined_w250_s20_GTR_F_I_G4_20250619_021313",
        "results/cleaned_alignment_combined_w250_s25_GTR_F_I_G4_20250619_204021",
    ]

    print("🧬 FOCUSED TREE GENOTYPE MAPPING")
    print("=" * 80)
    print("Target files: all_rooted_trees.*, rooted_trees_collection.*")
    print("Features: Keep accessions, remove duplicates, better file names")
    print()

    if not os.path.exists(fasta_file):
        print(f"❌ FASTA file not found: {fasta_file}")
        return

    all_processed_files = []

    for target_dir in target_dirs:
        results_dir = os.path.join(base_dir, target_dir)
        processed_files = process_results_directory(results_dir, fasta_file)
        if processed_files:
            all_processed_files.extend(processed_files)

    print(f"\n🎉 COMPLETE! Total files processed: {len(all_processed_files)}")
    print("\nProcessed files with improved genotype mappings:")
    for f in all_processed_files:
        print(f"   📄 {os.path.relpath(f, base_dir)}")


if __name__ == "__main__":
    main()
