#!/usr/bin/env python3
"""
Script to map genogroup names back to phylogenetic trees.
This script extracts genogroup information from FASTA headers and replaces
accession numbers in tree files with genogroup names.

Usage: python map_genogroups_to_trees.py
"""

import os
import re
import glob


def parse_fasta_headers(fasta_file):
    """
    Parse FASTA file to extract mapping between accession numbers and genogroup names.

    Args:
        fasta_file (str): Path to the FASTA file

    Returns:
        dictif __name__ == "__main__":
    # You can choose which function to run:
    # main()  # For processing individual window tree files
    # main_backbone_trees()  # For processing ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre
    # main_specific_tree_file()  # For processing a specific tree file in its own directory
    main_with_combined_sequences()  # For processing tree files with cleaned combined sequencesping of accession number to genogroup name
    """
    accession_to_genogroup = {}

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Example header: >KR074148.1 Norovirus GII/Hu/BRA/2004/GII.P7-GII.6/RS7842 RdRp and VP1 genes, partial cds [trimmed 15:15] [refined]
                header = line.strip()

                # Extract accession number (first part after >)
                accession_match = re.match(r">(\S+)", header)
                if not accession_match:
                    continue
                accession = accession_match.group(1)

                # Extract genogroup (pattern: GII.PX-GII.Y or GII.Pe-GII.Y or GII.Pg-GII.Y)
                genogroup_match = re.search(
                    r"(GII\.P[^-\s/]+(?:-GII\.\d+|17|12|15))", header
                )
                if genogroup_match:
                    genogroup = genogroup_match.group(1)
                    accession_to_genogroup[accession] = genogroup
                else:
                    print(f"Warning: Could not extract genogroup from header: {header}")

    return accession_to_genogroup


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


def process_tree_file(tree_file, accession_to_genogroup, output_dir):
    """
    Process a single tree file, replacing accession numbers with genogroup names.

    Args:
        tree_file (str): Path to the tree file
        accession_to_genogroup (dict): Mapping of accession to genogroup
        output_dir (str): Output directory for processed trees
    """
    with open(tree_file, "r") as f:
        tree_content = f.read().strip()

    # Replace accession numbers with genogroup names
    modified_tree = tree_content
    for accession, genogroup in accession_to_genogroup.items():
        # Replace accession with genogroup_accession format
        modified_tree = modified_tree.replace(accession, f"{genogroup}_{accession}")

    # Create output filename
    window_number = os.path.basename(os.path.dirname(tree_file))
    output_filename = f"window_{window_number}_genogroups.treefile"
    output_path = os.path.join(output_dir, output_filename)

    # Write modified tree
    with open(output_path, "w") as f:
        f.write(modified_tree)

    print(f"Processed tree for window {window_number} -> {output_filename}")

    return output_path


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

    print(f"Created mapping file: {mapping_file}")
    return mapping_file


def create_genogroup_only_trees(tree_files, accession_to_genogroup, output_dir):
    """
    Create trees with only genogroup names (no accession numbers).

    Args:
        tree_files (list): List of tree file paths
        accession_to_genogroup (dict): Mapping of accession to genogroup
        output_dir (str): Output directory
    """
    genogroup_only_dir = os.path.join(output_dir, "genogroup_only")
    os.makedirs(genogroup_only_dir, exist_ok=True)

    for tree_file in tree_files:
        with open(tree_file, "r") as f:
            tree_content = f.read().strip()

        # Replace accession numbers with genogroup names only
        modified_tree = tree_content
        for accession, genogroup in accession_to_genogroup.items():
            modified_tree = modified_tree.replace(accession, genogroup)

        # Create output filename
        window_number = os.path.basename(os.path.dirname(tree_file))
        output_filename = f"window_{window_number}_genogroups_only.treefile"
        output_path = os.path.join(genogroup_only_dir, output_filename)

        # Write modified tree
        with open(output_path, "w") as f:
            f.write(modified_tree)

        print(f"Created genogroup-only tree for window {window_number}")


def process_all_backbone_trees(tree_file, accession_to_genogroup, output_dir):
    """
    Process the ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre file, replacing accession numbers with genogroup names.

    Args:
        tree_file (str): Path to the ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre file
        accession_to_genogroup (dict): Mapping of accession to genogroup
        output_dir (str): Output directory for processed trees
    """
    with open(tree_file, "r") as f:
        lines = f.readlines()

    # Process each tree line
    modified_lines_with_accession = []
    modified_lines_genogroup_only = []

    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue

        # Replace accession numbers with genogroup_accession format
        modified_tree_with_accession = line
        for accession, genogroup in accession_to_genogroup.items():
            modified_tree_with_accession = modified_tree_with_accession.replace(
                accession, f"{genogroup}_{accession}"
            )
        modified_lines_with_accession.append(modified_tree_with_accession)

        # Replace accession numbers with genogroup names only
        modified_tree_genogroup_only = line
        for accession, genogroup in accession_to_genogroup.items():
            modified_tree_genogroup_only = modified_tree_genogroup_only.replace(
                accession, genogroup
            )
        modified_lines_genogroup_only.append(modified_tree_genogroup_only)

    # Create output filenames
    base_filename = os.path.splitext(os.path.basename(tree_file))[0]

    # Write modified tree with genogroup_accession format
    output_path_with_accession = os.path.join(
        output_dir, f"{base_filename}_with_genogroups.tre"
    )
    with open(output_path_with_accession, "w") as f:
        for line in modified_lines_with_accession:
            f.write(line + "\n")

    # Write modified tree with genogroup names only
    output_path_genogroup_only = os.path.join(
        output_dir, f"{base_filename}_genogroups_only.tre"
    )
    with open(output_path_genogroup_only, "w") as f:
        for line in modified_lines_genogroup_only:
            f.write(line + "\n")

    print("Processed ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre:")
    print(f"  - Created: {output_path_with_accession}")
    print(f"  - Created: {output_path_genogroup_only}")
    print(f"  - Processed {len(modified_lines_with_accession)} tree lines")

    return output_path_with_accession, output_path_genogroup_only


def main():
    """Main function to process all tree files."""

    # Configuration
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    fasta_file = os.path.join(base_dir, "data", "refined_alignment.fasta")
    results_dir = os.path.join(
        base_dir, "results_nf_backbone_midpoint_300_2", "tree_reconstruction", "iqtree"
    )
    output_dir = os.path.join(base_dir, "genogroup_mapped_trees")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Parse FASTA file to get accession to genogroup mapping
    print("Parsing FASTA file for genogroup information...")
    accession_to_genogroup = parse_fasta_headers(fasta_file)
    print(f"Found {len(accession_to_genogroup)} accession-genogroup mappings")

    # Print summary of genogroups found
    genogroups = set(accession_to_genogroup.values())
    print(f"Unique genogroups found: {sorted(genogroups)}")

    # Create mapping file
    create_mapping_file(accession_to_genogroup, output_dir)

    # Find all tree files
    tree_pattern = os.path.join(results_dir, "*", "*.treefile")
    tree_files = glob.glob(tree_pattern)
    print(f"Found {len(tree_files)} tree files to process")

    if not tree_files:
        print("No tree files found!")
        return

    # Process each tree file
    print("\nProcessing tree files...")
    processed_files = []
    for tree_file in sorted(tree_files):
        try:
            output_path = process_tree_file(
                tree_file, accession_to_genogroup, output_dir
            )
            processed_files.append(output_path)
        except Exception as e:
            print(f"Error processing {tree_file}: {e}")

    # Create genogroup-only trees
    print("\nCreating genogroup-only trees...")
    create_genogroup_only_trees(tree_files, accession_to_genogroup, output_dir)

    # Summary
    print("\nSummary:")
    print(f"- Processed {len(processed_files)} tree files")
    print(f"- Output directory: {output_dir}")
    print(f"- Mapping file created with {len(accession_to_genogroup)} entries")
    print(f"- Found {len(genogroups)} unique genogroups")

    # Show first few entries of mapping
    print("\nFirst 10 accession-genogroup mappings:")
    for i, (accession, genogroup) in enumerate(
        sorted(accession_to_genogroup.items())[:10]
    ):
        print(f"  {accession} -> {genogroup}")


def main_backbone_trees():
    """Main function to process the ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre file specifically."""

    # Configuration
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    fasta_file = os.path.join(base_dir, "data", "refined_alignment.fasta")
    backbone_tree_file = os.path.join(
        base_dir, "results", "ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre"
    )
    output_dir = os.path.join(base_dir, "genogroup_mapped_trees")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Check if files exist
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}")
        return

    if not os.path.exists(backbone_tree_file):
        print(f"Error: Backbone tree file not found: {backbone_tree_file}")
        return

    # Parse FASTA file to get accession to genogroup mapping
    print("Parsing FASTA file for genogroup information...")
    accession_to_genogroup = parse_fasta_headers(fasta_file)
    print(f"Found {len(accession_to_genogroup)} accession-genogroup mappings")

    # Print summary of genogroups found
    genogroups = set(accession_to_genogroup.values())
    print(f"Unique genogroups found: {sorted(genogroups)}")

    # Create mapping file
    create_mapping_file(accession_to_genogroup, output_dir)

    # Process the backbone tree file
    print("\nProcessing ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre file...")
    try:
        output_with_accession, output_genogroup_only = process_all_backbone_trees(
            backbone_tree_file, accession_to_genogroup, output_dir
        )

        # Summary
        print("\nSummary:")
        print(f"- Output directory: {output_dir}")
        print(f"- Mapping file created with {len(accession_to_genogroup)} entries")
        print(f"- Found {len(genogroups)} unique genogroups")
        print("- Created tree files:")
        print(
            f"  * With genogroups and accessions: {os.path.basename(output_with_accession)}"
        )
        print(f"  * Genogroups only: {os.path.basename(output_genogroup_only)}")

        # Show first few entries of mapping
        print("\nFirst 10 accession-genogroup mappings:")
        for i, (accession, genogroup) in enumerate(
            sorted(accession_to_genogroup.items())[:10]
        ):
            print(f"  {accession} -> {genogroup}")

    except Exception as e:
        print(f"Error processing backbone tree file: {e}")


def main_specific_tree_file():
    """Main function to process a specific tree file in its own directory."""

    # Configuration - specific tree file
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    fasta_file = os.path.join(base_dir, "data", "refined_alignment.fasta")
    tree_file = os.path.join(
        base_dir,
        "results",
        "ingroup_only_alignment_w300_s10_GTR_F_I_G4_20250618_214333",
        "all_rooted_trees.newick",
    )

    # Output directory - same as the tree file directory
    output_dir = os.path.dirname(tree_file)

    # Check if files exist
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}")
        return

    if not os.path.exists(tree_file):
        print(f"Error: Tree file not found: {tree_file}")
        return

    # Parse FASTA file to get accession to genogroup mapping
    print("Parsing FASTA file for genogroup information...")
    accession_to_genogroup = parse_fasta_headers(fasta_file)
    print(f"Found {len(accession_to_genogroup)} accession-genogroup mappings")

    # Print summary of genogroups found
    genogroups = set(accession_to_genogroup.values())
    print(f"Unique genogroups found: {sorted(genogroups)}")

    # Create mapping file in the same directory as the tree file
    create_mapping_file(accession_to_genogroup, output_dir)

    # Process the tree file
    print(f"\nProcessing {os.path.basename(tree_file)} file...")
    try:
        output_with_accession, output_genogroup_only = process_all_backbone_trees(
            tree_file, accession_to_genogroup, output_dir
        )

        # Summary
        print("\nSummary:")
        print(f"- Output directory: {output_dir}")
        print(f"- Mapping file created with {len(accession_to_genogroup)} entries")
        print(f"- Found {len(genogroups)} unique genogroups")
        print("- Created tree files:")
        print(
            f"  * With genogroups and accessions: {os.path.basename(output_with_accession)}"
        )
        print(f"  * Genogroups only: {os.path.basename(output_genogroup_only)}")

        # Show first few entries of mapping
        print("\nFirst 10 accession-genogroup mappings:")
        for i, (accession, genogroup) in enumerate(
            sorted(accession_to_genogroup.items())[:10]
        ):
            print(f"  {accession} -> {genogroup}")

    except Exception as e:
        print(f"Error processing tree file: {e}")


def load_enhanced_mappings(base_dir):
    """Load optimized mappings from the enhanced analysis."""
    mapping_file = os.path.join(
        base_dir, "enhanced_genotype_analysis", "phylogenetic_tree_mapping.tsv"
    )

    if not os.path.exists(mapping_file):
        print(f"Warning: Enhanced mapping file not found: {mapping_file}")
        print(
            "Run enhanced_norovirus_analysis.py first to generate optimized mappings."
        )
        return {}

    optimized_mapping = {}
    with open(mapping_file, "r") as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                accession = parts[0]
                optimized_name = parts[1]
                optimized_mapping[accession] = optimized_name

    return optimized_mapping


def main_with_combined_sequences():
    """Main function to process tree files using the cleaned combined sequences with enhanced mappings."""

    # Configuration - using cleaned combined sequences
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    fasta_file = os.path.join(base_dir, "data", "cleaned_alignment_combined.fasta")
    tree_file = os.path.join(
        base_dir,
        "results",
        "cleaned_alignment_combined_w200_s20_GTR_F_I_G4_20250619_073748",
        "all_rooted_trees.newick",
    )

    # Output directory - same as the tree file directory
    output_dir = os.path.dirname(tree_file)

    # Check if files exist
    if not os.path.exists(fasta_file):
        print(f"Error: Cleaned FASTA file not found: {fasta_file}")
        print(
            "Please run the cleaning script first to create the cleaned alignment file."
        )
        return

    if not os.path.exists(tree_file):
        print(f"Error: Tree file not found: {tree_file}")
        return

    # Load enhanced optimized mappings
    print("Loading enhanced mappings for optimized phylogenetic tree display...")
    optimized_mapping = load_enhanced_mappings(base_dir)

    # Parse FASTA file to get accession to genogroup mapping (fallback)
    print("Parsing cleaned alignment file for genogroup information...")
    accession_to_genogroup = parse_fasta_headers_enhanced(fasta_file)
    print(f"Found {len(accession_to_genogroup)} accession-genogroup mappings")
    print(f"Found {len(optimized_mapping)} optimized display names")

    # Print summary of genogroups found
    genogroups = set(accession_to_genogroup.values())
    print(f"Unique genogroups found: {sorted(genogroups)}")

    # Create mapping file in the same directory as the tree file
    create_mapping_file(accession_to_genogroup, output_dir)

    # Process the tree file with enhanced mappings
    print(f"\nProcessing {os.path.basename(tree_file)} file with enhanced mappings...")
    try:
        # Create enhanced tree files
        with open(tree_file, "r") as f:
            lines = f.readlines()

        # Process each tree line
        modified_lines_with_accession = []
        modified_lines_genogroup_only = []
        modified_lines_optimized = []

        for line in lines:
            line = line.strip()
            if not line:
                continue

            # Version 1: Replace accession numbers with genogroup_accession format
            modified_tree_with_accession = line
            for accession, genogroup in accession_to_genogroup.items():
                modified_tree_with_accession = modified_tree_with_accession.replace(
                    accession, f"{genogroup}_{accession}"
                )
            modified_lines_with_accession.append(modified_tree_with_accession)

            # Version 2: Replace accession numbers with genogroup names only
            modified_tree_genogroup_only = line
            for accession, genogroup in accession_to_genogroup.items():
                modified_tree_genogroup_only = modified_tree_genogroup_only.replace(
                    accession, genogroup
                )
            modified_lines_genogroup_only.append(modified_tree_genogroup_only)

            # Version 3: Replace with optimized display names
            modified_tree_optimized = line
            for accession in accession_to_genogroup.keys():
                if accession in optimized_mapping:
                    optimized_name = optimized_mapping[accession]
                    modified_tree_optimized = modified_tree_optimized.replace(
                        accession, optimized_name
                    )
                else:
                    # Fallback to genogroup if optimized name not available
                    genogroup = accession_to_genogroup[accession]
                    modified_tree_optimized = modified_tree_optimized.replace(
                        accession, genogroup
                    )
            modified_lines_optimized.append(modified_tree_optimized)

        # Create output filenames
        base_filename = os.path.splitext(os.path.basename(tree_file))[0]

        # Write enhanced tree files
        output_paths = {}

        # Standard genogroup_accession format
        output_with_accession = os.path.join(
            output_dir, f"{base_filename}_with_genogroups.tre"
        )
        with open(output_with_accession, "w") as f:
            for line in modified_lines_with_accession:
                f.write(line + "\n")
        output_paths["with_accession"] = output_with_accession

        # Genogroup names only
        output_genogroup_only = os.path.join(
            output_dir, f"{base_filename}_genogroups_only.tre"
        )
        with open(output_genogroup_only, "w") as f:
            for line in modified_lines_genogroup_only:
                f.write(line + "\n")
        output_paths["genogroup_only"] = output_genogroup_only

        # Optimized display names
        output_optimized = os.path.join(
            output_dir, f"{base_filename}_optimized_display.tre"
        )
        with open(output_optimized, "w") as f:
            for line in modified_lines_optimized:
                f.write(line + "\n")
        output_paths["optimized"] = output_optimized

        # Summary
        print("\nSummary:")
        print(f"- Output directory: {output_dir}")
        print(f"- Mapping file created with {len(accession_to_genogroup)} entries")
        print(f"- Found {len(genogroups)} unique genogroups")
        print("- Created tree files:")
        print(
            f"  * With genogroups and accessions: {os.path.basename(output_with_accession)}"
        )
        print(f"  * Genogroups only: {os.path.basename(output_genogroup_only)}")
        print(f"  * Optimized display names: {os.path.basename(output_optimized)}")

        # Show first few entries of mapping
        print("\nFirst 10 accession-genogroup mappings from cleaned dataset:")
        for i, (accession, genogroup) in enumerate(
            sorted(accession_to_genogroup.items())[:10]
        ):
            optimized_name = optimized_mapping.get(accession, "N/A")
            print(f"  {accession} -> {genogroup} -> {optimized_name}")

        # Show some paper sequences specifically (if any remain after cleaning)
        print("\nSample mappings from paper sequences (if present in cleaned data):")
        paper_sequences = [
            acc for acc in accession_to_genogroup.keys() if acc.startswith(("KY", "MF"))
        ]
        for i, acc in enumerate(sorted(paper_sequences)[:5]):
            optimized_name = optimized_mapping.get(acc, "N/A")
            print(f"  {acc} -> {accession_to_genogroup[acc]} -> {optimized_name}")

        return output_paths

    except Exception as e:
        print(f"Error processing tree file: {e}")
        return None


if __name__ == "__main__":
    # You can choose which function to run:
    # main()  # For processing individual window tree files
    # main_backbone_trees()  # For processing ALL_BACKBONE_MIDPOINT_TREES_SIMPLE.tre
    # main_specific_tree_file()  # For processing a specific tree file in its own directory
    main_with_combined_sequences()  # For processing tree files with combined sequences
