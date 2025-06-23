#!/usr/bin/env python3

import os
import sys
import subprocess
import json
from datetime import datetime


def log_message(message, log_file):
    """Log message with timestamp."""
    timestamp = datetime.now().isoformat()
    with open(log_file, "a") as f:
        f.write(f"[{timestamp}] {message}\n")
    print(f"[{timestamp}] {message}")


def clean_alignment_headers(alignment_file, output_file):
    """Clean alignment headers to ensure compatibility."""
    with open(alignment_file, "r") as f:
        lines = f.readlines()

    with open(output_file, "w") as f:
        for line in lines:
            if line.startswith(">"):
                # Extract accession number only
                import re

                match = re.search(r">([^|\s]+)", line)
                if match:
                    f.write(f">{match.group(1)}\n")
                else:
                    f.write(line)
            else:
                f.write(line)

    return output_file


def run_iqtree_backbone(alignment_file, output_prefix, model, cpus):
    """Run IQ-TREE to generate backbone tree from full alignment."""
    try:
        # Clean alignment
        clean_alignment = clean_alignment_headers(
            alignment_file, f"{output_prefix}_clean.fasta"
        )

        # Run IQ-TREE with automatic model selection
        cmd = [
            "iqtree2",
            "-s",
            clean_alignment,
            "--prefix",
            output_prefix,
            "-m",
            model,
            "--threads",
            str(cpus),
            "-B",
            "1000",  # Bootstrap support
            "--quiet",
        ]

        log_message(
            f"Running IQ-TREE for backbone tree: {' '.join(cmd)}", "backbone_tree.log"
        )

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            log_message(
                "IQ-TREE backbone analysis completed successfully", "backbone_tree.log"
            )
            return f"{output_prefix}.treefile"
        else:
            log_message(f"IQ-TREE failed: {result.stderr}", "backbone_tree.log")
            return None

    except Exception as e:
        log_message(f"Error running IQ-TREE: {str(e)}", "backbone_tree.log")
        return None


def apply_midpoint_rooting(tree_file, output_file):
    """Apply midpoint rooting to backbone tree."""
    try:
        import dendropy

        # Load tree
        tree = dendropy.Tree.get(path=tree_file, schema="newick")

        # Apply midpoint rooting
        tree.reroot_at_midpoint()

        # Write rooted tree
        tree.write(path=output_file, schema="newick")

        log_message("Midpoint rooting applied to backbone tree", "backbone_tree.log")
        return True

    except Exception as e:
        log_message(f"Error applying midpoint rooting: {str(e)}", "backbone_tree.log")
        return False


def main():
    # Parse command line arguments
    if len(sys.argv) != 4:
        print("Usage: generate_backbone_tree.py <alignment_file> <model> <cpus>")
        sys.exit(1)

    alignment_file = sys.argv[1]
    model = sys.argv[2] if sys.argv[2] != "null" else "MFP"
    cpus = int(sys.argv[3])

    log_file = "backbone_tree.log"

    # Initialize log
    log_message("=== Backbone Tree Generation ===", log_file)
    log_message(f"Input alignment: {alignment_file}", log_file)
    log_message(f"Model: {model}", log_file)

    success = False

    try:
        # Step 1: Generate backbone tree with IQ-TREE
        backbone_tree = run_iqtree_backbone(
            alignment_file, "backbone_tree", model, cpus
        )

        if backbone_tree and os.path.exists(backbone_tree):
            log_message("Backbone tree generation: SUCCESS", log_file)

            # Step 2: Apply midpoint rooting to backbone tree
            if apply_midpoint_rooting(backbone_tree, "backbone_tree_rooted.treefile"):
                log_message("Backbone midpoint rooting: SUCCESS", log_file)

                # Use rooted backbone as final backbone tree
                os.rename("backbone_tree_rooted.treefile", "backbone_tree.treefile")
                success = True
            else:
                log_message(
                    "Backbone midpoint rooting: FAILED - using unrooted tree", log_file
                )
                success = True  # Still proceed with unrooted backbone
        else:
            log_message("Backbone tree generation: FAILED", log_file)

    except Exception as e:
        log_message(f"Backbone analysis failed with exception: {str(e)}", log_file)

    # Create output file even if failed (for pipeline continuity)
    if not success:
        log_message("Creating empty backbone tree due to analysis failure", log_file)
        with open("backbone_tree.treefile", "w") as f:
            f.write("# Backbone tree generation failed\n")

    # Final status
    status = "SUCCESS" if success else "FAILED"
    log_message(f"=== BACKBONE TREE GENERATION: {status} ===", log_file)

    # Generate analysis summary
    summary = {
        "analysis_type": "backbone_tree_generation",
        "timestamp": datetime.now().isoformat(),
        "status": status,
        "parameters": {
            "model": model,
            "threads": cpus,
            "bootstrap": 1000,
            "rooting_method": "midpoint",
        },
    }

    with open(f"{log_file}.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Generate versions file
    with open("versions.yml", "w") as f:
        f.write('"GENERATE_BACKBONE_TREE":\n')

        # Get IQ-TREE version
        try:
            iqtree_version = subprocess.check_output(
                ["iqtree2", "--version"], stderr=subprocess.STDOUT, text=True
            )
            version_line = iqtree_version.split("\n")[0]
            version = (
                version_line.split()[2] if len(version_line.split()) >= 3 else "unknown"
            )
            f.write(f'    iqtree: "{version}"\n')
        except subprocess.CalledProcessError:
            f.write('    iqtree: "unknown"\n')

        # Get Python version
        f.write(f'    python: "{sys.version.split()[0]}"\n')

        # Check for dendropy
        try:
            import dendropy

            f.write(f'    dendropy: "{dendropy.__version__}"\n')
        except ImportError:
            f.write('    dendropy: "not available"\n')


if __name__ == "__main__":
    main()
