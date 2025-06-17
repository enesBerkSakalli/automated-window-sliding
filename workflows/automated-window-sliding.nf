/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import java.text.SimpleDateFormat
import java.util.Date

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def summary_params = paramsSummaryMap(workflow)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { SlidingWindow      } from '../modules/local/sliding_window'
include { ModelSelection     } from '../modules/local/model_finder'
include { IQTREE2            } from '../modules/local/tree_reconstruction/iqtree2'
include { RAXMLNG            } from '../modules/local/tree_reconstruction/raxmlng'
include { CollectTrees       } from '../modules/local/collect_trees'
include { MAD_ROOTING        } from '../modules/local/mad_rooting'
include { COLLECT_ROOTED_TREES } from '../modules/local/collect_rooted_trees'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AutomatedWindowSliding {

    // Construct unique output directory with model information
    def inputFile = file(params.input)
    def inputFileBasename = inputFile.getBaseName()
    def timestamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
    def modelName = params.model ? params.model.replaceAll('[+/]', '_') : 'auto'
    def unique_outdir = "${params.outdir}/${inputFileBasename}_w${params.window_size}_s${params.step_size}_${modelName}_${timestamp}"
    new File(unique_outdir).mkdirs()

    // Preserve original alignment in output directory for reproducibility
    process CopyOriginalAlignment {
        publishDir unique_outdir, mode: 'copy'
        
        input:
        path(alignment)
        
        output:
        path("original_alignment_${inputFileBasename}.fasta")
        
        script:
        """
        cp ${alignment} original_alignment_${inputFileBasename}.fasta
        """
    }
    
    // Generate analysis metadata for reproducibility
    process GenerateAnalysisMetadata {
        publishDir unique_outdir, mode: 'copy'
        
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
Pipeline version: ${workflow.manifest.version}
Command line: ${workflow.commandLine}
Work directory: ${workflow.workDir}
Launch directory: ${workflow.launchDir}

## System Information
Container: ${workflow.container}
Profile: ${workflow.profile}
EOF
        """
    }
    
    // Copy the original alignment to results
    input_alignment = Channel.fromPath(params.input)
    original_alignment_copy = CopyOriginalAlignment(input_alignment)
    analysis_metadata = GenerateAnalysisMetadata()

    // -- MAIN WORKFLOW --

    // check if a file with custom windows was provided -> if not the dummy file NO_FILE is used
    // needs to be done because there is no option for optional input channels
    if ( params.window_file ) {
        custom_windows_file = Channel.fromPath(params.window_file)
    } else {
        custom_windows_file = Channel.fromPath("${projectDir}/assets/NO_FILE")
    }

    // SlidingWindow call
    sliding_window_output = SlidingWindow(input_alignment, custom_windows_file)


    // collect a overview of the windows
    sliding_window_output.log.collectFile(storeDir: unique_outdir)
    
    // collect the removed sequences for each window
    sliding_window_output.removed_sequences.collectFile(storeDir: unique_outdir)

    // caluculate the relative length for each window of the total length of all windows for resource allocation
    window_lengths = sliding_window_output.log.splitCsv(sep: '\t', header: true).map({ row -> tuple(row.name, row.win_len.toFloat())})
    total_length = window_lengths.map({it[1]}).sum()
    window_lengths = window_lengths.combine(total_length).map({tuple(it[0], it[1]/it[2])})

    alignment_windows = sliding_window_output.alignments

    num_windows = alignment_windows.flatten().count()

    // assign the model channel
    if ( params.model != null && !params.model.isEmpty() ) {
        // create a channel from the model specified with --model and store it in a file
        evolutionary_model = Channel.from([params.model]).collectFile(name: 'model.txt')
    } else {
        if ( params.model_finder_splits ) {
            // add the relative length of the window to the input for model finder for resource allocation
            model_finder_input = alignment_windows.flatten().map({tuple(it.baseName, it)}).join(window_lengths).map({tuple(it[1], it[2])})
        } else {
            // use dummy channel to satisfy input cardinality
            model_finder_input = input_alignment.combine(Channel.from([1]))
        }
        // ModelSelection module call - assuming no publishDir changes needed within it for now
        evolutionary_model = ModelSelection(model_finder_input, num_windows)

        // output tsv file with models for each window
        evolutionary_model.model_tuple.map({"${it[0]}\t${it[1].text}"}).collectFile(name: 'models.txt', storeDir: unique_outdir, sort: {it.split('\t')[0].isInteger()?it.split('\t')[0].toInteger():it.split('\t')[0]})
        evolutionary_model.log.collectFile(storeDir: "${unique_outdir}/model_finder_logs/logs")
        evolutionary_model.iqtree.collectFile(storeDir: "${unique_outdir}/model_finder_logs/iqtree")
    }

    // prepare the input for tree reconstruction -> combine alignment windows with evolutionary model
    if ( params.model != null && !params.model.isEmpty() ) {
        tree_input = alignment_windows.flatten().combine(evolutionary_model.map({it.text.trim()}))
    } else {
        if ( params.model_finder_splits ) {
            alignment_windows = alignment_windows.flatten().map( { tuple(it.baseName, it) } )
            tree_input = alignment_windows.join(evolutionary_model.model_tuple).map( {tuple(it[1], it[2].text.trim())} )
        } else {
            tree_input = alignment_windows.flatten().combine(evolutionary_model.model_tuple.map( {it[1].text.trim()}))
        }
    }

    // add the relative length of the window to the input for model finder for resource allocation
    tree_input = tree_input.map({tuple(it[0].baseName, it[0], it[1])}).join(window_lengths).map({tuple(it[1], it[2], it[3])})

    // tree input now a tuple with 3 entries (alignment_file, model_str, relative_length)

    // choose tree reconstruction program that should be used
    // to add a new phylogenetic program define new process in modules/local and import into this script
    // add the new phylogenetic program as a new branch in this if-else block. Also add parameter as valid pattern in nextflow_schema.json
    // new process should have a output channel containing the reconstructed trees. Assign this output channel to a new channel called treefiles.

    if ( params.phylo_method == "iqtree2" ) {
        // IQTREE2 module call - will update later if module needs unique_outdir param
        tree = IQTREE2(tree_input, num_windows)
        treefiles = tree.treefile
        contrees = tree.consensus_tree
        tree.iqtree.collectFile(storeDir: "${unique_outdir}/tree_reconstruction_logs/iqtree_files")
    } else if ( params.phylo_method == "raxml-ng" ) {
        // RAXMLNG module call - will update later if module needs unique_outdir param
        tree = RAXMLNG(tree_input, num_windows)
        treefiles = tree.best_tree
        contrees = tree.consensus_tree
    } else {
        throw new Error("Phylogenetic Program'${params.phylo_method}' not supported")
    }

    tree.log.collectFile(storeDir: "${unique_outdir}/tree_reconstruction_logs/logs")
    treefiles = treefiles.collect().map( {tuple("best_trees", it)} )
    contrees = contrees.collect().map( {tuple("consensus_trees", it)} )
    trees = treefiles.concat(contrees)

    // MAD rooting for best trees (if enabled)
    if (params.mad_rooting) {
        // Create a map of alignment files by their base names
        alignment_map = alignment_windows.flatten()
            .map { alignment_file -> 
                def window_name = alignment_file.simpleName
                tuple(window_name, alignment_file)
            }
        
        // Extract only the best trees for rooting (not consensus trees)
        best_trees_for_rooting = tree.treefile
            .flatten()
            .map { tree_file -> 
                def window_name = tree_file.simpleName
                tuple(window_name, tree_file)
            }
        
        // Combine trees with their corresponding alignments
        trees_with_alignments = best_trees_for_rooting.join(alignment_map)
            .map { window_name, tree_file, alignment_file ->
                tuple(window_name, tree_file, alignment_file)
            }
        
        // Apply MAD rooting to best trees with alignments
        rooted_trees = MAD_ROOTING(trees_with_alignments, unique_outdir)
        
        // Collect rooted trees and generate reports
        rooted_collection = COLLECT_ROOTED_TREES(
            rooted_trees.rooted_trees.map{it[1]}.collect(),
            rooted_trees.logs.collect(),
            unique_outdir,
            params.output_format
        )
    }

    // sort the trees by the midpoint position of the window
    if (params.window_file) {
        // create hashmap to get midpoint position based on the filename
        midpoints = sliding_window_output.log.splitCsv(sep: '\t', header: true).map({tuple(it.name, it.mid)}).collect(flat: false).collectEntries(elem -> [elem[0], elem[1].toInteger()])
        trees = trees.map({tuple(it[0], it[1].toSorted({elem -> midpoints.get(elem.simpleName).value}))})
    } else {
        trees = trees.map({tuple(it[0], it[1].toSorted({elem -> elem.simpleName.toInteger()}))})
    }

    // CollectTrees module call - will update later if module needs unique_outdir param
    CollectTrees(trees)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
