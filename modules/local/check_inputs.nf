#!/usr/bin/env nextflow

// Input validation workflow
workflow CHECK_INPUTS {
    def required_files = [
        params.alignment_file
    ]
    
    required_files.each { file_path ->
        if (!file(file_path).exists()) {
            error "Required input file missing: ${file_path}"
        } else {
            log.info "✓ Found input file: ${file_path}"
        }
    }
}

workflow {
    CHECK_INPUTS()
}
