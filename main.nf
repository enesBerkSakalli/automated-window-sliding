#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    automated-window-sliding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ggruber193/automated-window-sliding
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { validateParameters; paramsHelp } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AutomatedWindowSliding } from './workflows/automated-window-sliding'

//
// WORKFLOW: Run main automated-window-sliding analysis pipeline
//
workflow AutomatedWindowSlidingWorkflow {
    AutomatedWindowSliding()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    // Print help message if needed
    if (params.help) {
        def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
        def String command = "nextflow run ${workflow.manifest.name} --input <ALIGNMENT> --outdir <OUTDIR> --window_size <int> --step_size <int>"
        log.info(
            logo + command + """
        
        NOTE: All results are automatically organized in a structured directory under the 'results/' folder.
        Each run creates a timestamped subdirectory for reproducibility and version control.
        
        """.stripIndent()
        )
        System.exit(0)
    }

    // Print workflow version and exit on --version
    if (params.version) {
        def String workflow_version = NfcoreTemplate.version(workflow)
        log.info("${workflow.manifest.name} ${workflow_version}")
        System.exit(0)
    }

    // Ensure the input file exists before processing
    if (!file(params.input).exists()) {
        error("Input alignment file does not exist: ${params.input}")
    }

    AutomatedWindowSlidingWorkflow()
}
