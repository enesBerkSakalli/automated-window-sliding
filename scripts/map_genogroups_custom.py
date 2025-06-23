#!/usr/bin/env python3
"""
Custom script to map genogroup names to phylogenetic trees in the cleaned alignment results.
This script processes the results from: results/cleaned_alignment_combined_w250_s10_GTR_F_I_G4_20250619_103929

Usage: python map_genogroups_custom.py
"""

import os
import re
import glob


def parse_fasta_headers_enhanced(fasta_file):
    """
    Enhanced version to parse FASTA file to extract mapping between accession numbers and genogroup names.
    Works with both original and new sequence header formats.

    Args:
        fasta_file (str): Path to the FASTA file

    Returns:
        dict: Mapping of accession number to genogroup name
    """
    accession_to_genogroup = {}

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()

                # Extract accession number (first part after >)
                accession_match = re.match(r">(\S+)", header)
                if not accession_match:
                    continue
                accession = accession_match.group(1)

                # Try multiple patterns for genogroup extraction
                genogroup = None

                # Pattern 1: GII.PX-GII.Y format (most common)
                genogroup_match = re.search(
                    r"(GI{1,2}\.P[^-\s/]+(?:-GI{1,2}\.\d+|17|12|15))", header
                )
                if genogroup_match:
                    genogroup = genogroup_match.group(1)

                # Pattern 2: Simple GII.X format
                elif re.search(r"(GI{1,2}\.\d+)", header):
                    simple_match = re.search(r"(GI{1,2}\.\d+)", header)
                    genogroup = simple_match.group(1)

                # Pattern 3: GII.Pe, GII.Pg formats
                elif re.search(r"(GI{1,2}\.P[a-z]+)", header):
                    pe_match = re.search(r"(GI{1,2}\.P[a-z]+)", header)
                    genogroup = pe_match.group(1)

                if genogroup:
                    accession_to_genogroup[accession] = genogroup
                else:
                    print(f"Warning: Could not extract genogroup from header: {header}")

    return accession_to_genogroup


def process_single_tree_file(tree_file, accession_to_genogroup, output_file_suffix=""):
    """
    Process a single tree file, replacing accession numbers with genogroup names.

    Args:
        tree_file (str): Path to the tree file
        accession_to_genogroup (dict): Mapping of accession to genogroup
        output_file_suffix (str): Suffix to add to output filename

    Returns:
        tuple: (path_with_genogroups, path_genogroups_only)
    """
    with open(tree_file, "r") as f:
        tree_content = f.read().strip()

    # Replace accession numbers with genogroup_accession format
    modified_tree_with_accession = tree_content
    for accession, genogroup in accession_to_genogroup.items():
        modified_tree_with_accession = modified_tree_with_accession.replace(
            accession, f"{genogroup}_{accession}"
        )

    # Replace accession numbers with genogroup names only
    modified_tree_genogroup_only = tree_content
    for accession, genogroup in accession_to_genogroup.items():
        modified_tree_genogroup_only = modified_tree_genogroup_only.replace(
            accession, genogroup
        )

    # Create output filenames
    base_dir = os.path.dirname(tree_file)
    base_name = os.path.splitext(os.path.basename(tree_file))[0]

    output_with_accession = os.path.join(
        base_dir, f"{base_name}_with_genogroups{output_file_suffix}.treefile"
    )
    output_genogroup_only = os.path.join(
        base_dir, f"{base_name}_genogroups_only{output_file_suffix}.treefile"
    )

    # Write modified trees
    with open(output_with_accession, "w") as f:
        f.write(modified_tree_with_accession)

    with open(output_genogroup_only, "w") as f:
        f.write(modified_tree_genogroup_only)

    return output_with_accession, output_genogroup_only


def create_mapping_file(accession_to_genogroup, output_dir):
    """
    Create a mapping file showing the relationship between accession numbers and genogroups.

    Args:
        accession_to_genogroup (dict): Mapping of accession to genogroup
        output_dir (str): Output directory
    """
    mapping_file = os.path.join(output_dir, "accession_genogroup_mapping.txt")

    with open(mapping_file, "w") as f:
        f.write("Accession_Number\tGenogroup\n")
        for accession, genogroup in sorted(accession_to_genogroup.items()):
            f.write(f"{accession}\t{genogroup}\n")

    print(f"✅ Created mapping file: {mapping_file}")
    return mapping_file


def process_cleaned_alignment_results():
    """
    Main function to process tree files from the cleaned alignment results.
    Specifically targets: results/cleaned_alignment_combined_w250_s10_GTR_F_I_G4_20250619_103929
    """

    # Configuration
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    fasta_file = os.path.join(base_dir, "data", "cleaned_alignment_combined.fasta")
    results_dir = os.path.join(
        base_dir,
        "results",
        "cleaned_alignment_combined_w250_s10_GTR_F_I_G4_20250619_103929",
    )

    # Check if files exist
    if not os.path.exists(fasta_file):
        print(f"❌ Error: Cleaned FASTA file not found: {fasta_file}")
        return

    if not os.path.exists(results_dir):
        print(f"❌ Error: Results directory not found: {results_dir}")
        return

    print("🔬 GENOGROUP MAPPING FOR CLEANED ALIGNMENT RESULTS")
    print("=" * 80)
    print(f"📁 Results Directory: {os.path.basename(results_dir)}")
    print()

    # Parse FASTA file to get accession to genogroup mapping
    print("📋 Parsing cleaned alignment file for genogroup information...")
    accession_to_genogroup = parse_fasta_headers_enhanced(fasta_file)
    print(f"✅ Found {len(accession_to_genogroup)} accession-genogroup mappings")

    # Print summary of genogroups found
    genogroups = set(accession_to_genogroup.values())
    print(f"🧬 Unique genogroups found: {len(genogroups)}")
    print(f"   {sorted(genogroups)}")
    print()

    # Create mapping file in the results directory
    create_mapping_file(accession_to_genogroup, results_dir)
    print()

    # Process different types of tree files
    processed_files = {}

    # 1. Process main tree files in root directory
    print("🌳 Processing main tree files...")
    main_tree_files = [
        "all_rooted_trees.newick",
        "all_constraint_trees.newick",
        "rooted_trees_collection.newick",
    ]

    for tree_filename in main_tree_files:
        tree_path = os.path.join(results_dir, tree_filename)
        if os.path.exists(tree_path):
            try:
                with_geo, geo_only = process_single_tree_file(
                    tree_path, accession_to_genogroup
                )
                processed_files[tree_filename] = {
                    "with_genogroups": with_geo,
                    "genogroups_only": geo_only,
                }
                print(f"   ✅ Processed: {tree_filename}")
            except Exception as e:
                print(f"   ❌ Error processing {tree_filename}: {e}")
        else:
            print(f"   ⚠️  File not found: {tree_filename}")

    print()

    # 2. Process individual window trees in backbone_midpoint directory
    print("🪟 Processing individual window trees...")
    backbone_dir = os.path.join(results_dir, "backbone_midpoint")

    if os.path.exists(backbone_dir):
        # Find all rooted tree files
        rooted_tree_pattern = os.path.join(backbone_dir, "*_rooted.treefile")
        rooted_tree_files = glob.glob(rooted_tree_pattern)

        # Find all constraint tree files
        constraint_tree_pattern = os.path.join(backbone_dir, "*_constraint.treefile")
        constraint_tree_files = glob.glob(constraint_tree_pattern)

        print(f"   📊 Found {len(rooted_tree_files)} rooted tree files")
        print(f"   📊 Found {len(constraint_tree_files)} constraint tree files")

        # Process rooted trees
        processed_rooted = 0
        for tree_file in sorted(rooted_tree_files):
            try:
                window_num = os.path.basename(tree_file).split("_")[0]
                with_geo, geo_only = process_single_tree_file(
                    tree_file, accession_to_genogroup
                )
                processed_rooted += 1
                if processed_rooted <= 3:  # Show first few
                    print(f"   ✅ Processed rooted tree for window {window_num}")
            except Exception as e:
                print(f"   ❌ Error processing {os.path.basename(tree_file)}: {e}")

        if processed_rooted > 3:
            print(f"   ✅ ... and {processed_rooted - 3} more rooted trees")

        # Process constraint trees
        processed_constraint = 0
        for tree_file in sorted(constraint_tree_files):
            try:
                window_num = os.path.basename(tree_file).split("_")[0]
                with_geo, geo_only = process_single_tree_file(
                    tree_file, accession_to_genogroup
                )
                processed_constraint += 1
                if processed_constraint <= 3:  # Show first few
                    print(f"   ✅ Processed constraint tree for window {window_num}")
            except Exception as e:
                print(f"   ❌ Error processing {os.path.basename(tree_file)}: {e}")

        if processed_constraint > 3:
            print(f"   ✅ ... and {processed_constraint - 3} more constraint trees")

        processed_files["individual_windows"] = {
            "rooted_trees": processed_rooted,
            "constraint_trees": processed_constraint,
        }
    else:
        print(f"   ⚠️  Backbone directory not found: {backbone_dir}")

    print()

    # Summary
    print("📈 PROCESSING SUMMARY")
    print("-" * 50)
    print(f"✅ Results directory: {os.path.basename(results_dir)}")
    print(f"✅ Sequences processed: {len(accession_to_genogroup)}")
    print(f"✅ Unique genogroups: {len(genogroups)}")
    print(
        f"✅ Main tree files processed: {len([f for f in processed_files if f != 'individual_windows'])}"
    )

    if "individual_windows" in processed_files:
        print(
            f"✅ Individual rooted trees: {processed_files['individual_windows']['rooted_trees']}"
        )
        print(
            f"✅ Individual constraint trees: {processed_files['individual_windows']['constraint_trees']}"
        )

    print()
    print("📁 OUTPUT FILES CREATED:")
    print("   • accession_genogroup_mapping.txt - Complete mapping reference")
    print("   • *_with_genogroups.treefile - Trees with genogroup_accession format")
    print("   • *_genogroups_only.treefile - Trees with genogroup names only")
    print()

    # Show sample mappings
    print("🔍 SAMPLE ACCESSION-GENOGROUP MAPPINGS:")
    print("-" * 50)
    for i, (accession, genogroup) in enumerate(
        sorted(accession_to_genogroup.items())[:10]
    ):
        print(f"   {accession} → {genogroup}")
    if len(accession_to_genogroup) > 10:
        print(f"   ... and {len(accession_to_genogroup) - 10} more mappings")

    print()
    print("🎯 NEXT STEPS:")
    print("-" * 50)
    print("   📊 Use *_genogroups_only.treefile for visualization")
    print("   🔬 Use *_with_genogroups.treefile for detailed analysis")
    print("   📋 Reference accession_genogroup_mapping.txt for full details")

    print()
    print("=" * 80)
    print("🏆 GENOGROUP MAPPING COMPLETED SUCCESSFULLY!")
    print("All tree files have been processed with genogroup annotations.")
    print("=" * 80)


if __name__ == "__main__":
    process_cleaned_alignment_results()
