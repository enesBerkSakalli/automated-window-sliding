#!/usr/bin/env python3
"""
Comprehensive Norovirus Table Generator and Tree Mapper
Detective-style analysis to extract genotypes and create publication-ready tables and tree mappings.

This script:
1. Extracts genotype data from FASTA files with enhanced pattern matching
2. Validates against known publications
3. Creates comprehensive LaTeX tables with publication integration
4. Applies optimized genotype names to phylogenetic trees
5. Handles complex genotype formats and edge cases

Author: Enhanced AI Detective System
Date: 2025
"""

import os
import re
import glob
from collections import defaultdict, Counter
from pathlib import Path
import json


class NorovirusDetective:
    """Enhanced detective class for comprehensive norovirus genotype analysis."""

    def __init__(self):
        self.setup_patterns()
        self.setup_publications()
        self.setup_standardizations()

    def setup_patterns(self):
        """Setup enhanced regex patterns for genotype detection."""
        self.genotype_patterns = [
            # Pattern 1: Standard dual-typing with dash (GII.PX-GII.Y)
            (r"(GI{1,2}\.P\w+(?:\d+)?[-]GI{1,2}\.\d+(?:_\w+)?)", "dual_dash", 10),
            # Pattern 2: Dual-typing with underscore (GII.PX_GII.Y) - for MF sequences
            (r"(GI{1,2}\.P\w+(?:\d+)?[_]GI{1,2}\.\d+(?:_\w+)?)", "dual_underscore", 9),
            # Pattern 3: Standard dual-typing with tilde (GII.PX~GII.Y)
            (r"(GI{1,2}\.P\w+(?:\d+)?[~]GI{1,2}\.\d+(?:_\w+)?)", "dual_tilde", 8),
            # Pattern 4: MF sequence special format (GII.Pe_GII.4_Sydney)
            (r"(GI{1,2}\.P[a-z]+[_]GI{1,2}\.\d+[_]\w+)", "mf_special", 7),
            # Pattern 5: Direct in description "Norovirus GII.X"
            (r"Norovirus\s+(GI{1,2}\.\d+(?:\w+)?)", "norovirus_prefix", 6),
            # Pattern 6: P-types with letters
            (r"(GI{1,2}\.P[a-z]+(?:\d+)?)", "p_type_letters", 5),
            # Pattern 7: Variant names with spaces
            (r"(GI{1,2}\.\d+\s+\w+(?:\s+\w+)?)", "variant_names", 4),
            # Pattern 8: Standard genotypes
            (r"(GI{1,2}\.\d+(?:\w+)?)", "standard_genotypes", 3),
            # Pattern 9: P-types only
            (r"(GI{1,2}\.P\d+)", "p_type_only", 2),
        ]

    def setup_publications(self):
        """Setup publication database for cross-referencing."""
        self.publications = {
            "barreira_2017": {
                "citation": "Barreira et~al. (2017)",
                "full_citation": "Barreira, D.M.P.G., et al. Detection and molecular characterization of the novel recombinant norovirus GII.P16-GII.4 Sydney in southeastern Brazil in 2016. PLoS ONE 12(12): e0189504.",
                "doi": "10.1371/journal.pone.0189504",
                "genotypes": [
                    "GII.P16-GII.4",
                    "GII.Pe-GII.4",
                    "GII.P17-GII.17",
                    "GII.P16-GII.3",
                    "GII.Pg-GII.1",
                ],
                "years": [2015, 2016],
                "location": "Brazil",
            },
            "current_study": {
                "citation": "Current study",
                "full_citation": "Current phylogenetic analysis of norovirus strains from Brazil (2004-2016)",
                "genotypes": ["various"],
                "years": list(range(2004, 2017)),
                "location": "Brazil",
            },
            "fumian_2016": {
                "citation": "Fumian et~al. (2016)",
                "full_citation": "Fumian, T.M., et al. Norovirus Recombinant Strains Isolated from Gastroenteritis Outbreaks in Southern Brazil, 2004-2011. PLoS One. 2016 Apr 26;11(4):e0145391.",
                "doi": "10.1371/journal.pone.0145391",
                "genotypes": [
                    "GII.P7-GII.6",
                    "GII.P21-GII.3",
                    "GII.P4-GII.4",
                    "GII.P7-GII.14",
                    "GII.P31-GII.17",
                    "GII.P13-GII.17",
                    "GII.P21-GII.21",
                    "GII.P2-GII.2",
                    "GII.P12-GII.12",
                    "GII.P21-GII.13",
                    "GII.P16-GII.3",
                    "GII.P15-GII.15",
                ],
                "years": list(range(2004, 2012)),
                "location": "Brazil",
            },
        }

    def setup_standardizations(self):
        """Setup genotype standardization mappings."""
        self.genotype_standardization = {
            "GII.Pe": "GII.P31",  # Updated nomenclature
            "GII.Pg": "GII.P12",
            "GII.Pa": "GII.P4",
            "GII.PNA": "GII.PNA1",
        }

        # Phylogenetic tree display optimization
        self.display_optimization = {
            "max_length": 25,
            "abbreviations": {
                "Sydney": "Syd",
                "New Orleans": "NO",
                "Norovirus": "NV",
                "partial": "part",
                "complete": "comp",
                "genome": "gen",
                "strain": "st",
                "nonstructural": "nonstruct",
                "polyprotein": "polyprot",
            },
        }

    def extract_comprehensive_info(self, header):
        """
        Enhanced comprehensive genotype extraction with detective-level analysis.

        Args:
            header (str): FASTA header line

        Returns:
            dict: Comprehensive genotype information
        """
        info = {
            "accession": None,
            "dual_type": None,
            "genotype": None,
            "p_type": None,
            "variant": None,
            "year": None,
            "country": None,
            "strain_id": None,
            "display_name": None,
            "publication": None,
            "confidence": 0,
        }

        # Extract accession number
        acc_match = re.match(r">(\S+)", header)
        if acc_match:
            info["accession"] = acc_match.group(1)

        # Extract year
        year_patterns = [r"/(\d{4})/", r"_(\d{4})", r"\b(\d{4})\b"]
        for pattern in year_patterns:
            match = re.search(pattern, header)
            if match:
                year = int(match.group(1))
                if 1970 <= year <= 2030:
                    info["year"] = year
                    break

        # Extract country
        country_match = re.search(r"/([A-Z]{2,3})/", header)
        if country_match:
            info["country"] = country_match.group(1)

        # Extract strain ID
        strain_patterns = [
            r"/([A-Z0-9_]+)\s+(?:RdRp|nonstructural|genes)",
            r"/([A-Z0-9_]+)$",
            r"([A-Z0-9_]{4,})",
        ]
        for pattern in strain_patterns:
            match = re.search(pattern, header)
            if match:
                info["strain_id"] = match.group(1)
                break

        # Enhanced genotype extraction with confidence scoring
        best_match = None
        best_confidence = 0

        for pattern, pattern_type, confidence in self.genotype_patterns:
            match = re.search(pattern, header)
            if match and confidence > best_confidence:
                extracted = match.group(1)
                best_match = extracted
                best_confidence = confidence
                info["confidence"] = confidence

                # Handle dual-typing formats
                if "_" in extracted or "-" in extracted or "~" in extracted:
                    extracted = extracted.replace("_", "-").replace("~", "-")
                    info["dual_type"] = extracted
                    parts = re.split(r"[-]", extracted)
                    if len(parts) == 2:
                        info["p_type"] = parts[0]
                        info["genotype"] = parts[1]

                # Handle P-types
                elif extracted.startswith(("GI.P", "GII.P")):
                    info["p_type"] = extracted
                    genotype_match = re.search(r"GI{1,2}\.\d+(?!\s*P)", header)
                    if genotype_match:
                        info["genotype"] = genotype_match.group(0)
                        info["dual_type"] = f"{extracted}-{info['genotype']}"

                # Handle direct genotypes
                else:
                    info["genotype"] = extracted

        # Additional MF sequence handling
        if not info["dual_type"] and "strain" in header.lower():
            mf_pattern = re.search(r"(GI{1,2}\.P[a-z]+)[_-](GI{1,2}\.\d+)", header)
            if mf_pattern:
                info["p_type"] = mf_pattern.group(1)
                info["genotype"] = mf_pattern.group(2)
                info["dual_type"] = f"{info['p_type']}-{info['genotype']}"
                info["confidence"] = 8

        # Standardize P-types
        if info["p_type"] and info["p_type"] in self.genotype_standardization:
            old_p_type = info["p_type"]
            info["p_type"] = self.genotype_standardization[old_p_type]
            if info["dual_type"]:
                info["dual_type"] = info["dual_type"].replace(
                    old_p_type, info["p_type"]
                )

        # Assign publication based on genotype and year
        info["publication"] = self.assign_publication(info)

        # Create optimized display name
        info["display_name"] = self.create_optimized_display_name(info)

        return info

    def assign_publication(self, info):
        """Assign publication based on genotype patterns and years."""
        if not info.get("year"):
            return "Current study"

        year = info["year"]
        dual_type = info.get("dual_type", "") or ""  # Handle None case

        # Barreira et al. 2017 patterns
        if year in [2015, 2016]:
            barreira_genotypes = [
                "GII.P16-GII.4",
                "GII.Pe-GII.4",
                "GII.P17-GII.17",
                "GII.P16-GII.3",
                "GII.Pg-GII.1",
            ]
            if any(gt in dual_type for gt in barreira_genotypes):
                return "Barreira et~al. (2017)"

        # Fumian et al. 2016 patterns (KR074148-KR074189 accession numbers)
        if 2004 <= year <= 2011:
            fumian_genotypes = [
                "GII.P7-GII.6",
                "GII.P21-GII.3",
                "GII.P4-GII.4",
                "GII.P7-GII.14",
                "GII.P31-GII.17",
                "GII.P13-GII.17",
                "GII.P21-GII.21",
                "GII.P2-GII.2",
                "GII.P12-GII.12",
                "GII.P21-GII.13",
                "GII.P16-GII.3",
                "GII.P15-GII.15",
            ]
            if any(gt in dual_type for gt in fumian_genotypes):
                return "Fumian et~al. (2016)"

        return "Current study"

    def create_optimized_display_name(self, info):
        """Create optimized phylogenetic tree display name."""
        components = []

        # Add genotype information
        if info["dual_type"]:
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
                info["strain_id"][:8]
                if len(info["strain_id"]) > 8
                else info["strain_id"]
            )
            components.append(strain_short)

        # Create final name with length limit
        display_name = "_".join(components)
        if len(display_name) > self.display_optimization["max_length"]:
            display_name = display_name[: self.display_optimization["max_length"]]

        return display_name

    def analyze_fasta_file(self, fasta_file):
        """Analyze FASTA file and extract comprehensive genotype information."""
        print(f"Detective Analysis: Investigating {fasta_file}")

        sequences_data = []

        with open(fasta_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith(">"):
                    header = line.strip()
                    info = self.extract_comprehensive_info(header)
                    info["line_number"] = line_num
                    sequences_data.append(info)

        print(f"Detective Report: Analyzed {len(sequences_data)} sequences")
        return sequences_data

    def create_comprehensive_table(self, sequences_data, output_dir):
        """Create comprehensive LaTeX table with publication integration."""

        # Group data by year and genotype
        grouped_data = defaultdict(lambda: defaultdict(list))

        for seq in sequences_data:
            if seq["year"] and seq["dual_type"]:
                year = seq["year"]
                genotype = seq["dual_type"]
                publication = seq["publication"]
                grouped_data[year][(genotype, publication)].append(seq)

        # Create LaTeX table content
        latex_content = self.generate_latex_table(grouped_data)

        # Write table file
        table_file = os.path.join(output_dir, "comprehensive_norovirus_table.tex")
        with open(table_file, "w") as f:
            f.write(latex_content)

        print(f"Comprehensive table created: {table_file}")
        return table_file

    def generate_latex_table(self, grouped_data):
        """Generate comprehensive LaTeX table content."""

        latex_content = """\\begin{table}[htbp]
\t\\centering
\t\\caption{Distribution of norovirus sequenced strains by year, genotype, frequency, and source (comprehensive analysis).}
\t\\scriptsize
\t\\begin{tabular}{p{0.6cm}p{3.2cm}p{1.8cm}p{2.5cm}}
\t\t\\toprule
\t\t\\textbf{Year} & \\textbf{Norovirus region ORF1-ORF2} & \\textbf{Strains (\\%) N=XX} & \\textbf{Publication} \\\\
\t\t\\midrule
"""

        total_strains = sum(
            len(seqs)
            for year_data in grouped_data.values()
            for seqs in year_data.values()
        )

        # Sort years
        for year in sorted(grouped_data.keys()):
            year_data = grouped_data[year]

            # Sort genotypes by frequency (descending)
            sorted_genotypes = sorted(
                year_data.items(), key=lambda x: len(x[1]), reverse=True
            )

            for (genotype, publication), sequences in sorted_genotypes:
                count = len(sequences)
                percentage = (count / total_strains) * 100

                # Format genotype for LaTeX (replace - with ~)
                genotype_latex = genotype.replace("-", "~")

                latex_content += f"\t\t{year:<13} & {genotype_latex:<45} & {count} ({percentage:.1f})  & {publication:<30} \\\\\n"

        latex_content += """\t\t\\bottomrule
\t\\end{tabular}
\t\\label{tab:comprehensive_norovirus_strains}
\\end{table}"""

        # Update total count in caption
        latex_content = latex_content.replace("N=XX", f"N={total_strains}")

        return latex_content

    def apply_to_tree_files(self, sequences_data, results_dirs):
        """Apply optimized genotype names to phylogenetic tree files."""

        # Create accession to display name mapping
        accession_mapping = {}
        for seq in sequences_data:
            if seq["accession"] and seq["display_name"]:
                accession_mapping[seq["accession"]] = seq["display_name"]

        print(f"Detective Mapping: Created {len(accession_mapping)} accession mappings")

        # Apply to each results directory
        for results_dir in results_dirs:
            if os.path.exists(results_dir):
                print(f"Applying mappings to: {results_dir}")
                self.process_tree_directory(results_dir, accession_mapping)
            else:
                print(f"Warning: Directory not found: {results_dir}")

    def process_tree_directory(self, results_dir, accession_mapping):
        """Process all tree files in a results directory."""

        # Find all tree files
        tree_patterns = ["*.treefile", "*.newick", "*.nexus", "*.tre"]
        tree_files = []

        for pattern in tree_patterns:
            tree_files.extend(
                glob.glob(os.path.join(results_dir, "**", pattern), recursive=True)
            )

        print(f"Found {len(tree_files)} tree files in {results_dir}")

        # Process each tree file
        for tree_file in tree_files:
            self.update_tree_file(tree_file, accession_mapping)

    def update_tree_file(self, tree_file, accession_mapping):
        """Update a single tree file with genotype mappings."""

        try:
            # Read original tree file
            with open(tree_file, "r") as f:
                tree_content = f.read()

            # Track replacements
            replacements_made = 0

            # Apply mappings
            for accession, display_name in accession_mapping.items():
                if accession in tree_content:
                    tree_content = tree_content.replace(accession, display_name)
                    replacements_made += 1

            # Write updated tree file if changes were made
            if replacements_made > 0:
                output_file = tree_file.replace(".", "_genotype_mapped.")
                with open(output_file, "w") as f:
                    f.write(tree_content)
                print(
                    f"Updated {tree_file} -> {output_file} ({replacements_made} replacements)"
                )

        except Exception as e:
            print(f"Error processing {tree_file}: {e}")

    def generate_summary_report(self, sequences_data, output_dir):
        """Generate comprehensive summary report."""

        report_content = []
        report_content.append("COMPREHENSIVE NOROVIRUS DETECTIVE ANALYSIS REPORT")
        report_content.append("=" * 60)
        report_content.append("")

        # Basic statistics
        total_sequences = len(sequences_data)
        sequences_with_genotype = len([s for s in sequences_data if s["genotype"]])
        sequences_with_dual_type = len([s for s in sequences_data if s["dual_type"]])

        report_content.append(f"Total sequences analyzed: {total_sequences}")
        report_content.append(
            f"Sequences with genotype: {sequences_with_genotype} ({sequences_with_genotype / total_sequences * 100:.1f}%)"
        )
        report_content.append(
            f"Sequences with dual-type: {sequences_with_dual_type} ({sequences_with_dual_type / total_sequences * 100:.1f}%)"
        )
        report_content.append("")

        # Year distribution
        years = [s["year"] for s in sequences_data if s["year"]]
        if years:
            report_content.append(f"Year range: {min(years)} - {max(years)}")
            year_counts = Counter(years)
            report_content.append("Year distribution:")
            for year in sorted(year_counts.keys()):
                count = year_counts[year]
                report_content.append(f"  {year}: {count} sequences")
        report_content.append("")

        # Genotype distribution
        genotypes = [s["dual_type"] for s in sequences_data if s["dual_type"]]
        if genotypes:
            genotype_counts = Counter(genotypes)
            report_content.append("Top genotypes:")
            for genotype, count in genotype_counts.most_common(10):
                percentage = (count / len(genotypes)) * 100
                report_content.append(f"  {genotype}: {count} ({percentage:.1f}%)")
        report_content.append("")

        # Publication distribution
        publications = [s["publication"] for s in sequences_data if s["publication"]]
        if publications:
            pub_counts = Counter(publications)
            report_content.append("Publication distribution:")
            for pub, count in pub_counts.most_common():
                percentage = (count / len(publications)) * 100
                report_content.append(f"  {pub}: {count} ({percentage:.1f}%)")
        report_content.append("")

        # Confidence analysis
        confidences = [s["confidence"] for s in sequences_data if s["confidence"] > 0]
        if confidences:
            avg_confidence = sum(confidences) / len(confidences)
            report_content.append(
                f"Average extraction confidence: {avg_confidence:.1f}"
            )

        # Write report
        report_file = os.path.join(output_dir, "detective_analysis_report.txt")
        with open(report_file, "w") as f:
            f.write("\\n".join(report_content))

        print(f"Detective report generated: {report_file}")
        return report_file


def main():
    """Main execution function."""

    print("🔍 NOROVIRUS DETECTIVE SYSTEM ACTIVATED")
    print("=" * 50)

    # Initialize detective
    detective = NorovirusDetective()

    # Define paths
    fasta_file = "/Users/berksakalli/Projects/automated-window-sliding/data/cleaned_alignment_combined.fasta"
    output_dir = "/Users/berksakalli/Projects/automated-window-sliding/enhanced_genotype_analysis"

    # Results directories to process
    results_dirs = [
        "/Users/berksakalli/Projects/automated-window-sliding/results/cleaned_alignment_combined_w250_s10_GTR_F_I_G4_20250619_103929",
        "/Users/berksakalli/Projects/automated-window-sliding/results/cleaned_alignment_combined_w250_s20_GTR_F_I_G4_20250619_021313",
        "/Users/berksakalli/Projects/automated-window-sliding/results/cleaned_alignment_combined_w250_s25_GTR_F_I_G4_20250619_204021",
    ]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Analyze FASTA file
    print("\\n🔍 Step 1: Analyzing FASTA file...")
    sequences_data = detective.analyze_fasta_file(fasta_file)

    # Step 2: Create comprehensive table
    print("\\n📊 Step 2: Creating comprehensive table...")
    table_file = detective.create_comprehensive_table(sequences_data, output_dir)

    # Step 3: Apply to tree files
    print("\\n🌳 Step 3: Applying mappings to tree files...")
    detective.apply_to_tree_files(sequences_data, results_dirs)

    # Step 4: Generate summary report
    print("\\n📋 Step 4: Generating summary report...")
    report_file = detective.generate_summary_report(sequences_data, output_dir)

    print("\\n✅ DETECTIVE ANALYSIS COMPLETE!")
    print("Files generated:")
    print(f"  - Comprehensive table: {table_file}")
    print(f"  - Summary report: {report_file}")
    print(f"  - Updated tree files in: {len(results_dirs)} directories")


if __name__ == "__main__":
    main()
