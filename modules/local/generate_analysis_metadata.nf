process GenerateAnalysisMetadata {
    publishDir { unique_outdir }, mode: 'copy'
    input:
    val unique_outdir
    val inputFileBasename
    
    output:
    path("analysis_metadata.txt")
    
    script:
    """
    cat > analysis_metadata.txt << EOF
# Sliding Window Phylogenetic Analysis Metadata
# Generated: \$(date)

## Input Data
Original alignment file: ${params.input}
Original filename: ${inputFileBasename}.fasta

## Analysis Parameters
Window size: ${params.window_size}
Step size: ${params.step_size}
Evolutionary model: ${params.model ?: 'Auto-detected'}
MAD rooting enabled: ${params.mad_rooting}
Output format: ${params.output_format}

## Pipeline Information
Nextflow version: ${workflow.nextflow.version}
Pipeline version: ${workflow.manifest.version ?: 'unknown'}
Command line: ${workflow.commandLine}
Work directory: ${workflow.workDir}
Launch directory: ${workflow.launchDir}

## System Information
Container: ${workflow.container ?: 'none'}
Profile: ${workflow.profile}
EOF
    """
}
