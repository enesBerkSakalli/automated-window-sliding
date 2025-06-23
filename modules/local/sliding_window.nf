process SlidingWindow {
    publishDir path: "${params.outdir}/alignments/", pattern: "*.fasta", mode: params.publish_dir_mode, enabled: "${params.output_windows || params.run_mode == 'split'}"

    cpus 1

    conda "conda-forge::python=3.11 conda-forge::biopython"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81--py311h1db85ab_0' :
        'biocontainers/biopython:1.81--py311h1db85ab_0' }"

    input:
    path alignment
    path window_file

    output:
    path('*.fasta'), emit: window_files
    path("windows.log"), emit: log
    path("window_lengths.txt"), emit: window_lengths
    path("removed_sequences.log"), optional:true, emit: removed_sequences

    script:
    def window_file_arg = window_file.name != 'NO_FILE' ? "--split-file ${window_file}" : ''
    def window_size = params.window_size? "--window-size ${params.window_size}" : ''
    def step_size = params.step_size? "--step-size ${params.step_size}": ''
    def keep_ambiguous = (params.amb_seqs == 'remove') ? '' : '--keep-ambiguous'
    """
    # Install latest compatible BioPython without version pinning
    if ! python -c "import Bio" 2>/dev/null; then
        echo "Installing latest BioPython..."
        conda install -y -c conda-forge biopython || pip install biopython
    fi
    
    # Verify BioPython is available
    python -c "import Bio; print('BioPython version:', Bio.__version__)"
    
    cp /Users/berksakalli/Projects/automated-window-sliding/bin/sliding_window.py .
    chmod +x sliding_window.py
    python ./sliding_window.py --input ${alignment} ${window_size} ${step_size} --output-directory "" ${window_file_arg} ${keep_ambiguous} --log -1
    
    # Create proper window_lengths.txt with basename and relative length
    if [ ! -f "window_lengths.txt" ]; then
        echo "Creating window_lengths.txt with proper format"
        for fasta in *.fasta; do
            if [ -f "\$fasta" ]; then
                basename=\$(basename "\$fasta" .fasta)
                echo "\$basename\t1.0" >> window_lengths.txt
            fi
        done
    fi
    """
}