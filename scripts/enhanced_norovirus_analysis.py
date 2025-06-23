#!/usr/bin/env python3
"""
Enhanced Norovirus Genotype Analysis Script

This script analyzes norovirus sequences from FASTA files to:
1. Extract genotype information accurately following official nomenclature
2. Create frequency tables for LaTeX output
3. Optimize sequence names for phylogenetic tree display
4. Generate comprehensive mapping files

Based on norovirus classification system from:
- Chhabra et al. (2019) Updated classification of norovirus genogroups and genotypes
- Official norovirus typing tool nomenclature

Author: Enhanced for automated window sliding analysis
Date: 2025-06-20
"""

import os
import re
from collections import defaultdict, Counter
from datetime import datetime


class NorovirusGenotyper:
    """Class for comprehensive norovirus genotype analysis."""

    def __init__(self):
        """Initialize the genotyper with classification patterns."""
        # Official genotype patterns based on current nomenclature
        self.genotype_patterns = {
            # Standard dual-typing format: GII.PX-GII.Y
            "dual_typing": r"(GI{1,2}\.P\w+(?:\d+)?[-~]GI{1,2}\.\d+(?:\w+)?)",
            # Special P-types with letters
            "p_type_letters": r"(GI{1,2}\.P[a-z]+(?:\d+)?)",
            # Standard genotypes
            "genotype_only": r"(GI{1,2}\.\d+(?:\w+)?)",
            # Variants like "Sydney", "New Orleans"
            "variant_names": r"(GI{1,2}\.\d+\s+\w+(?:\s+\w+)?)",
            # P-types only
            "p_type_only": r"(GI{1,2}\.P\d+)",
        }

        # Mapping for standardizing genotype names
        self.genotype_standardization = {
            "GII.Pe": "GII.P31",  # According to updated nomenclature
            "GII.Pg": "GII.P12",
            "GII.Pa": "GII.P4",
            "GII.PNA": "GII.PNA1",
        }

        # Phylogenetic tree display optimization
        self.display_optimization = {
            "max_length": 20,  # Maximum length for tree labels
            "abbreviations": {
                "Sydney": "Syd",
                "New Orleans": "NO",
                "Norovirus": "NV",
                "partial": "part",
                "complete": "comp",
                "genome": "gen",
            },
        }

    def extract_year_from_header(self, header):
        """Extract year from FASTA header."""
        # Look for 4-digit year in various formats
        year_patterns = [
            r"/(\d{4})/",  # /2004/
            r"_(\d{4})",  # _2004
            r"\b(\d{4})\b",  # standalone 4-digit number
        ]

        for pattern in year_patterns:
            match = re.search(pattern, header)
            if match:
                year = int(match.group(1))
                # Validate year range (norovirus sequences typically from 1970s onward)
                if 1970 <= year <= 2030:
                    return year
        return None

    def extract_genotype_info(self, header):
        """
        Extract comprehensive genotype information from FASTA header.

        Returns:
            dict: Contains 'p_type', 'genotype', 'variant', 'dual_type', 'year'
        """
        info = {
            "p_type": None,
            "genotype": None,
            "variant": None,
            "dual_type": None,
            "year": None,
            "accession": None,
            "country": None,
            "strain_id": None,
        }

        # Extract accession number
        acc_match = re.match(r">(\S+)", header)
        if acc_match:
            info["accession"] = acc_match.group(1)

        # Extract year
        info["year"] = self.extract_year_from_header(header)

        # Extract country code
        country_match = re.search(r"/([A-Z]{2,3})/", header)
        if country_match:
            info["country"] = country_match.group(1)

        # Extract strain ID (last part of the path)
        strain_match = re.search(r"/([A-Z0-9_]+)\s+(?:RdRp|nonstructural)", header)
        if strain_match:
            info["strain_id"] = strain_match.group(1)

        # Enhanced genotype extraction patterns
        genotype_patterns = [
            # Pattern 1: Standard dual-typing with dash (GII.PX-GII.Y)
            r"(GI{1,2}\.P\w+(?:\d+)?[-]GI{1,2}\.\d+(?:_\w+)?)",
            # Pattern 2: Dual-typing with underscore (GII.PX_GII.Y) - for MF sequences
            r"(GI{1,2}\.P\w+(?:\d+)?[_]GI{1,2}\.\d+(?:_\w+)?)",
            # Pattern 3: Standard dual-typing with tilde (GII.PX~GII.Y)
            r"(GI{1,2}\.P\w+(?:\d+)?[~]GI{1,2}\.\d+(?:_\w+)?)",
            # Pattern 4: MF sequence special format (GII.Pe_GII.4_Sydney)
            r"(GI{1,2}\.P[a-z]+[_]GI{1,2}\.\d+[_]\w+)",
            # Pattern 5: Direct in description "Norovirus GII.X"
            r"Norovirus\s+(GI{1,2}\.\d+(?:\w+)?)",
            # Pattern 6: P-types with letters
            r"(GI{1,2}\.P[a-z]+(?:\d+)?)",
            # Pattern 7: Variant names
            r"(GI{1,2}\.\d+\s+\w+(?:\s+\w+)?)",
            # Pattern 8: Standard genotypes
            r"(GI{1,2}\.\d+(?:\w+)?)",
            # Pattern 9: P-types only
            r"(GI{1,2}\.P\d+)",
        ]

        # Try each pattern in order of specificity
        for pattern in genotype_patterns:
            match = re.search(pattern, header)
            if match:
                extracted = match.group(1)

                # Handle dual-typing formats
                if "_" in extracted or "-" in extracted or "~" in extracted:
                    # Convert underscore to dash for standardization
                    extracted = extracted.replace("_", "-").replace("~", "-")
                    info["dual_type"] = extracted
                    # Split into P-type and genotype
                    parts = re.split(r"[-]", extracted)
                    if len(parts) == 2:
                        info["p_type"] = parts[0]
                        info["genotype"] = parts[1]
                    break

                # Handle P-types
                elif extracted.startswith(("GI.P", "GII.P")):
                    info["p_type"] = extracted
                    # Look for corresponding genotype in the header
                    genotype_match = re.search(r"GI{1,2}\.\d+(?!\s*P)", header)
                    if genotype_match:
                        info["genotype"] = genotype_match.group(0)
                        info["dual_type"] = f"{extracted}-{info['genotype']}"

                # Handle direct genotypes
                else:
                    info["genotype"] = extracted
                    break

        # Additional parsing for complex headers like MF sequences
        if not info["dual_type"] and "strain" in header.lower():
            # Look for patterns like "GII.P16_GII.3" in strain name
            strain_pattern = re.search(r"(GI{1,2}\.P\w+[_-]GI{1,2}\.\d+)", header)
            if strain_pattern:
                dual_type = strain_pattern.group(1).replace("_", "-")
                info["dual_type"] = dual_type
                parts = dual_type.split("-")
                if len(parts) == 2:
                    info["p_type"] = parts[0]
                    info["genotype"] = parts[1]

        # Special handling for MF sequences with complex P-types like "GII.Pe"
        if not info["genotype"] or not info["p_type"]:
            # Handle MF sequences with format like "GII.Pe_GII.4_Sydney"
            mf_pattern = re.search(r"(GI{1,2}\.P[a-z]+)[_-](GI{1,2}\.\d+)", header)
            if mf_pattern:
                info["p_type"] = mf_pattern.group(1)
                info["genotype"] = mf_pattern.group(2)
                info["dual_type"] = f"{info['p_type']}-{info['genotype']}"

            # If still no genotype, try to extract from "Norovirus GII.X" at the beginning
            elif "Norovirus" in header:
                start_pattern = re.search(r"Norovirus\s+(GI{1,2}\.\d+)", header)
                if start_pattern:
                    info["genotype"] = start_pattern.group(1)

        # Standardize P-types according to updated nomenclature
        if info["p_type"] and info["p_type"] in self.genotype_standardization:
            info["p_type"] = self.genotype_standardization[info["p_type"]]

        return info

    def create_optimized_display_name(self, info):
        """
        Create optimized name for phylogenetic tree display.

        Args:
            info (dict): Genotype information from extract_genotype_info

        Returns:
            str: Optimized display name
        """
        components = []

        # Add genotype information
        if info["dual_type"]:
            # Use dual typing format but optimize
            dual_optimized = info["dual_type"]
            for full, abbrev in self.display_optimization["abbreviations"].items():
                dual_optimized = dual_optimized.replace(full, abbrev)
            components.append(dual_optimized)
        elif info["genotype"]:
            genotype_opt = info["genotype"]
            for full, abbrev in self.display_optimization["abbreviations"].items():
                genotype_opt = genotype_opt.replace(full, abbrev)
            components.append(genotype_opt)

        # Add year if available
        if info["year"]:
            components.append(str(info["year"]))

        # Add country if available
        if info["country"]:
            components.append(info["country"])

        # Add strain ID (shortened)
        if info["strain_id"]:
            strain_short = (
                info["strain_id"][:6]
                if len(info["strain_id"]) > 6
                else info["strain_id"]
            )
            components.append(strain_short)

        # Create final name
        display_name = "_".join(components)

        # Truncate if too long
        if len(display_name) > self.display_optimization["max_length"]:
            display_name = display_name[: self.display_optimization["max_length"]]

        return display_name

    def analyze_fasta_file(self, fasta_file):
        """
        Comprehensive analysis of FASTA file.

        Returns:
            dict: Analysis results including mappings and statistics
        """
        results = {
            "accession_to_info": {},
            "accession_to_display": {},
            "year_genotype_counts": defaultdict(lambda: defaultdict(int)),
            "genotype_counts": defaultdict(int),
            "year_counts": defaultdict(int),
            "errors": [],
        }

        with open(fasta_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith(">"):
                    header = line.strip()

                    try:
                        # Extract comprehensive information
                        info = self.extract_genotype_info(header)

                        if info["accession"]:
                            # Store full information
                            results["accession_to_info"][info["accession"]] = info

                            # Create optimized display name
                            display_name = self.create_optimized_display_name(info)
                            results["accession_to_display"][info["accession"]] = (
                                display_name
                            )

                            # Update statistics
                            if info["year"] and info["dual_type"]:
                                results["year_genotype_counts"][info["year"]][
                                    info["dual_type"]
                                ] += 1
                            elif info["year"] and info["genotype"]:
                                results["year_genotype_counts"][info["year"]][
                                    info["genotype"]
                                ] += 1

                            if info["dual_type"]:
                                results["genotype_counts"][info["dual_type"]] += 1
                            elif info["genotype"]:
                                results["genotype_counts"][info["genotype"]] += 1

                            if info["year"]:
                                results["year_counts"][info["year"]] += 1
                        else:
                            results["errors"].append(
                                f"Line {line_num}: Could not extract accession from {header}"
                            )

                    except Exception as e:
                        results["errors"].append(
                            f"Line {line_num}: Error processing header: {e}"
                        )

        return results

    def create_latex_table(self, analysis_results, existing_data=None):
        """
        Create extended LaTeX table with analysis results.

        Args:
            analysis_results (dict): Results from analyze_fasta_file
            existing_data (list): Existing table data to combine with

        Returns:
            str: LaTeX table code
        """
        # Prepare data for table
        table_data = []

        # Add existing data first (from the provided table)
        if existing_data:
            table_data.extend(existing_data)

        # Add new data from FASTA analysis
        year_genotype = analysis_results["year_genotype_counts"]
        total_sequences = sum(analysis_results["year_counts"].values())

        # Group by year and genotype
        for year in sorted(year_genotype.keys()):
            for genotype, count in sorted(year_genotype[year].items()):
                percentage = (count / total_sequences) * 100

                # Format genotype for display
                genotype_display = genotype.replace(
                    "-", "~"
                )  # Use tilde as in original table

                table_data.append(
                    {
                        "year": year,
                        "genotype": genotype_display,
                        "count": count,
                        "percentage": percentage,
                        "publication": "Current study",
                    }
                )

        # Create LaTeX table
        latex_table = (
            r"""
\begin{table}[htbp]
    \centering
    \caption{Distribution of norovirus sequenced strains by year, genotype, frequency, and source (supplementary data).}
    \scriptsize
    \begin{tabular}{p{0.6cm}p{3.2cm}p{1.8cm}p{2.5cm}}
        \toprule
        \textbf{Year} & \textbf{Norovirus region ORF1-ORF2} & \textbf{Strains (\%) N="""
            + str(total_sequences + 37)
            + r"""} & \textbf{Publication} \\
        \midrule"""
        )

        # Add existing data
        existing_entries = [
            (2015, "GII.P16~GII.3", 1, 2.7, "Barreira et~al. (2017)"),
            (2015, "GII.Pg~GII.1", 1, 2.7, "Barreira et~al. (2017)"),
            (2015, "GII.Pa~GII.4 Sydney 2012", 9, 24.3, "Barreira et~al. (2017)"),
            (2015, "GII.PNA~GII.4", 1, 2.7, "Barreira et~al. (2017)"),
            (2016, "GII.P17~GII.17", 3, 8.1, "Barreira et~al. (2017)"),
            (2016, "GII.P16~GII.4 Sydney", 22, 59.5, "Barreira et~al. (2017)"),
        ]

        # Combine and sort all data
        all_entries = []
        for entry in existing_entries:
            all_entries.append(entry)

        for data in table_data:
            all_entries.append(
                (
                    data["year"],
                    data["genotype"],
                    data["count"],
                    data["percentage"],
                    data["publication"],
                )
            )

        # Sort by year and genotype
        all_entries.sort(key=lambda x: (x[0], x[1]))

        # Add table rows
        for year, genotype, count, percentage, publication in all_entries:
            latex_table += f"""
        {year:<10} & {genotype:<35} & {count} ({percentage:.1f})  & {publication} \\\\"""

        latex_table += r"""
        \bottomrule
    \end{tabular}
    \label{tab:sup_norovirus_strains_extended}
\end{table}"""

        return latex_table

    def generate_summary_report(self, analysis_results, output_dir):
        """Generate comprehensive summary report."""
        report_file = os.path.join(output_dir, "norovirus_analysis_summary.txt")

        with open(report_file, "w") as f:
            f.write("NOROVIRUS GENOTYPE ANALYSIS SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            f.write(
                f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
            )

            # Basic statistics
            f.write("BASIC STATISTICS:\n")
            f.write("-" * 20 + "\n")
            f.write(
                f"Total sequences analyzed: {len(analysis_results['accession_to_info'])}\n"
            )
            f.write(
                f"Sequences with year information: {len([info for info in analysis_results['accession_to_info'].values() if info['year']])}\n"
            )
            f.write(
                f"Sequences with genotype information: {len([info for info in analysis_results['accession_to_info'].values() if info['dual_type'] or info['genotype']])}\n"
            )
            f.write(
                f"Unique genotypes found: {len(analysis_results['genotype_counts'])}\n"
            )
            f.write(
                f"Year range: {min(analysis_results['year_counts'].keys())} - {max(analysis_results['year_counts'].keys())}\n\n"
            )

            # Genotype distribution
            f.write("GENOTYPE DISTRIBUTION:\n")
            f.write("-" * 25 + "\n")
            for genotype, count in sorted(
                analysis_results["genotype_counts"].items(),
                key=lambda x: x[1],
                reverse=True,
            ):
                percentage = (count / len(analysis_results["accession_to_info"])) * 100
                f.write(f"{genotype:<25} {count:>3} ({percentage:>5.1f}%)\n")
            f.write("\n")

            # Year distribution
            f.write("YEAR DISTRIBUTION:\n")
            f.write("-" * 20 + "\n")
            for year in sorted(analysis_results["year_counts"].keys()):
                count = analysis_results["year_counts"][year]
                percentage = (count / len(analysis_results["accession_to_info"])) * 100
                f.write(f"{year}: {count:>3} sequences ({percentage:>5.1f}%)\n")
            f.write("\n")

            # Errors
            if analysis_results["errors"]:
                f.write("ERRORS AND WARNINGS:\n")
                f.write("-" * 25 + "\n")
                for error in analysis_results["errors"]:
                    f.write(f"{error}\n")

        print(f"Summary report created: {report_file}")
        return report_file

    def fix_missing_genotypes(self, analysis_results):
        """
        Fix missing genotype information by reconstructing from available data.

        Args:
            analysis_results (dict): Results from analyze_fasta_file

        Returns:
            dict: Updated analysis results with fixed genotypes
        """
        fixed_count = 0

        for accession, info in analysis_results["accession_to_info"].items():
            # If main genotype is missing but we have P-type and display name
            if (
                not info.get("dual_type")
                and not info.get("genotype")
                and info.get("p_type")
                and analysis_results["accession_to_display"].get(accession)
            ):
                display_name = analysis_results["accession_to_display"][accession]
                p_type = info["p_type"]

                # Extract genotype from display name
                # Pattern: GII.X_YEAR_COUNTRY
                genotype_match = re.match(r"(GI{1,2}\.\w+)_", display_name)
                if genotype_match:
                    inferred_genotype = genotype_match.group(1)

                    # Reconstruct dual-typing format
                    if p_type and inferred_genotype:
                        reconstructed_dual = f"{p_type}-{inferred_genotype}"
                        info["dual_type"] = reconstructed_dual
                        info["genotype"] = inferred_genotype

                        # Update statistics
                        if info["year"]:
                            analysis_results["year_genotype_counts"][info["year"]][
                                reconstructed_dual
                            ] += 1
                        analysis_results["genotype_counts"][reconstructed_dual] += 1

                        fixed_count += 1
                        print(f"Fixed genotype for {accession}: {reconstructed_dual}")

                # If no P-type but have inferred genotype from display name
                elif genotype_match and not p_type:
                    inferred_genotype = genotype_match.group(1)
                    info["genotype"] = inferred_genotype

                    # Update statistics
                    if info["year"]:
                        analysis_results["year_genotype_counts"][info["year"]][
                            inferred_genotype
                        ] += 1
                    analysis_results["genotype_counts"][inferred_genotype] += 1

                    fixed_count += 1
                    print(f"Added genotype for {accession}: {inferred_genotype}")

        print(f"Fixed {fixed_count} missing genotype entries")
        return analysis_results

    def create_comprehensive_mapping_table(self, analysis_results, output_dir):
        """
        Create a comprehensive mapping table with all available information.
        """
        comprehensive_file = os.path.join(
            output_dir, "comprehensive_genotype_mapping.tsv"
        )

        with open(comprehensive_file, "w") as f:
            f.write(
                "Accession\tDual_Type\tGenotype\tP_Type\tVariant\tYear\tCountry\tStrain_ID\tDisplay_Name\tStatus\n"
            )

            for acc, info in sorted(analysis_results["accession_to_info"].items()):
                display_name = analysis_results["accession_to_display"].get(acc, "")

                # Determine status
                status = "Complete"
                if not info.get("dual_type") and not info.get("genotype"):
                    status = "Missing_Genotype"
                elif not info.get("dual_type") and info.get("genotype"):
                    status = "Partial_P_Type_Missing"
                elif info.get("dual_type") and not info.get("p_type"):
                    status = "Partial_P_Type_Inferred"

                f.write(
                    f"{acc}\t{info.get('dual_type', '')}\t{info.get('genotype', '')}\t"
                    f"{info.get('p_type', '')}\t{info.get('variant', '')}\t"
                    f"{info.get('year', '')}\t{info.get('country', '')}\t"
                    f"{info.get('strain_id', '')}\t{display_name}\t{status}\n"
                )

        print(f"Created comprehensive mapping: {comprehensive_file}")
        return comprehensive_file


def main():
    """Main function to run the enhanced norovirus analysis."""

    # Configuration
    base_dir = "/Users/berksakalli/Projects/automated-window-sliding"
    fasta_file = os.path.join(base_dir, "data", "cleaned_alignment_combined.fasta")
    output_dir = os.path.join(base_dir, "enhanced_genotype_analysis")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Initialize genotyper
    genotyper = NorovirusGenotyper()

    # Check if FASTA file exists
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}")
        return

    print("Running enhanced norovirus genotype analysis...")
    print(f"Input file: {fasta_file}")
    print(f"Output directory: {output_dir}")
    print()

    # Analyze FASTA file
    print("Analyzing FASTA file...")
    analysis_results = genotyper.analyze_fasta_file(fasta_file)

    # Fix missing genotype information
    print("Fixing missing genotype information...")
    analysis_results = genotyper.fix_missing_genotypes(analysis_results)

    # Generate summary report
    print("Generating summary report...")
    genotyper.generate_summary_report(analysis_results, output_dir)

    # Create mapping files
    print("Creating mapping files...")

    # Standard accession to genotype mapping
    mapping_file = os.path.join(output_dir, "accession_genotype_mapping.tsv")
    with open(mapping_file, "w") as f:
        f.write("Accession\tGenotype\tP_Type\tYear\tCountry\tStrain_ID\tDisplay_Name\n")
        for acc, info in sorted(analysis_results["accession_to_info"].items()):
            display_name = analysis_results["accession_to_display"].get(acc, "")
            f.write(
                f"{acc}\t{info.get('dual_type', info.get('genotype', ''))}\t{info.get('p_type', '')}\t{info.get('year', '')}\t{info.get('country', '')}\t{info.get('strain_id', '')}\t{display_name}\n"
            )

    # Comprehensive mapping with status information
    genotyper.create_comprehensive_mapping_table(analysis_results, output_dir)

    # Phylogenetic tree mapping (accession -> optimized display name)
    tree_mapping_file = os.path.join(output_dir, "phylogenetic_tree_mapping.tsv")
    with open(tree_mapping_file, "w") as f:
        f.write("Original_Accession\tOptimized_Display_Name\n")
        for acc, display_name in sorted(
            analysis_results["accession_to_display"].items()
        ):
            f.write(f"{acc}\t{display_name}\n")

    # Create LaTeX table
    print("Creating extended LaTeX table...")
    latex_table = genotyper.create_latex_table(analysis_results)

    latex_file = os.path.join(output_dir, "extended_norovirus_table.tex")
    with open(latex_file, "w") as f:
        f.write(latex_table)

    print(f"\nAnalysis complete!")
    print(f"Files created:")
    print(
        f"  - Summary report: {os.path.basename(genotyper.generate_summary_report(analysis_results, output_dir))}"
    )
    print(f"  - Accession mapping: {os.path.basename(mapping_file)}")
    print(f"  - Tree mapping: {os.path.basename(tree_mapping_file)}")
    print(f"  - LaTeX table: {os.path.basename(latex_file)}")

    # Print quick summary
    print(f"\nQuick Summary:")
    print(f"  - Total sequences: {len(analysis_results['accession_to_info'])}")
    print(f"  - Unique genotypes: {len(analysis_results['genotype_counts'])}")
    print(
        f"  - Year range: {min(analysis_results['year_counts'].keys())} - {max(analysis_results['year_counts'].keys())}"
    )
    print(f"  - Top genotypes:")
    for genotype, count in sorted(
        analysis_results["genotype_counts"].items(), key=lambda x: x[1], reverse=True
    )[:5]:
        percentage = (count / len(analysis_results["accession_to_info"])) * 100
        print(f"    {genotype}: {count} ({percentage:.1f}%)")

    # Fix missing genotypes
    print("Fixing missing genotypes...")
    genotyper.fix_missing_genotypes(analysis_results)

    # Create comprehensive mapping table
    print("Creating comprehensive mapping table...")
    genotyper.create_comprehensive_mapping_table(analysis_results, output_dir)


if __name__ == "__main__":
    main()
