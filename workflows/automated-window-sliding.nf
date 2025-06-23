nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Note: Using fully-qualified names instead of import statements for Nextflow DSL2 compatibility

// include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'
// def summary_params = paramsSummaryMap(workflow)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { SlidingWindow } from '../modules/local/sliding_window'
include { ModelSelection } from '../modules/local/model_finder'
include { IQTREE2 } from '../modules/local/tree_reconstruction/iqtree2'
include { RAXMLNG } from '../modules/local/tree_reconstruction/raxmlng'
include { CollectTrees } from '../modules/local/collect_trees'
include { MAD_ROOTING } from '../modules/local/mad_rooting'
include { COLLECT_ROOTED_TREES } from '../modules/local/collect_rooted_trees'
include { GENERATE_BACKBONE_TREE } from '../modules/local/generate_backbone_tree'
include { BACKBONE_MIDPOINT_ROOTING } from '../modules/local/backbone_midpoint_rooting'
include { COLLECT_BACKBONE_MIDPOINT_TREES } from '../modules/local/collect_backbone_midpoint_trees'
include { GenerateAnalysisMetadata } from '../modules/local/generate_analysis_metadata'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AutomatedWindowSliding {

    // Define variables within workflow block (DSL2 requirement)
    def inputFile = file(params.input)
    def inputFileBasename = inputFile.getBaseName()
    def timestamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new java.util.Date())
    def modelName = params.model ? params.model.replaceAll('[+/]', '_') : 'auto'
    def main_results_dir = "results"
    def unique_outdir = "${main_results_dir}/${inputFileBasename}_w${params.window_size}_s${params.step_size}_${modelName}_${timestamp}"
    new File(unique_outdir).mkdirs()

    // Print results directory information
    log.info(
        """
    ===============================================================================
    AUTOMATED WINDOW SLIDING ANALYSIS
    ===============================================================================
    Input alignment: ${inputFile}
    Output directory: ${unique_outdir}
    Window size: ${params.window_size} bp
    Step size: ${params.step_size} bp
    Phylogenetic method: ${params.phylo_method}
    Model selection: ${params.model ?: 'automatic (' + params.model_criterion + ')'}
    Rooting: ${params.mad_rooting ? 'MAD' : 'none'}
    ===============================================================================
    """.stripIndent()
    )

    // Generate metadata
    GenerateAnalysisMetadata(unique_outdir, inputFileBasename)

    // Sliding window process - rename to avoid variable conflict
    def window_file_path = params.window_file ? file(params.window_file) : file('NO_FILE')
    SlidingWindow(inputFile, window_file_path)
    def alignment_windows = SlidingWindow.out.window_files
    def num_windows = alignment_windows.count()

    // Convert window_lengths file to channel for join operation
    def window_lengths = SlidingWindow.out.window_lengths
        .splitText()
        .map { line ->
            def parts = line.trim().split('\t')
            tuple(parts[0], parts[1].toFloat())
        }

    // Model selection
    def evolutionary_model
    if (params.model == null || params.model.isEmpty()) {
        def model_finder_input
        if (params.model_finder_splits) {
            model_finder_input = alignment_windows.flatten().map { tuple(it.baseName, it) }.join(window_lengths).map { tuple(it[1], it[2]) }
        }
        else {
            // use dummy channel to satisfy input cardinality - fix combine issue
            model_finder_input = Channel.of(inputFile).combine(Channel.of(1))
        }
        // ModelSelection module call - assuming no publishDir changes needed within it for now
        evolutionary_model = ModelSelection(model_finder_input, num_windows)

        // output tsv file with models for each window
        evolutionary_model.model_tuple.map { "${it[0]}\t${it[1].text}" }.collectFile(name: 'models.txt', storeDir: unique_outdir, sort: { it.split('\t')[0].isInteger() ? it.split('\t')[0].toInteger() : it.split('\t')[0] })
        evolutionary_model.log.collectFile(storeDir: "${unique_outdir}/model_finder_logs/logs")
        evolutionary_model.iqtree.collectFile(storeDir: "${unique_outdir}/model_finder_logs/iqtree")
    }
    else {
        evolutionary_model = Channel.of(params.model)
    }

    // prepare the input for tree reconstruction -> combine alignment windows with evolutionary model
    def tree_input
    if (params.model != null && !params.model.isEmpty()) {
        tree_input = alignment_windows.flatten().combine(evolutionary_model)
    }
    else {
        if (params.model_finder_splits) {
            def alignment_windows_tuple = alignment_windows.flatten().map { tuple(it.baseName, it) }
            tree_input = alignment_windows_tuple.join(evolutionary_model.model_tuple).map { tuple(it[1], it[2].text.trim()) }
        }
        else {
            tree_input = alignment_windows.flatten().combine(evolutionary_model.model_tuple.map { it[1].text.trim() })
        }
    }

    // add the relative length of the window to the input for model finder for resource allocation
    tree_input = tree_input.map { tuple(it[0].baseName, it[0], it[1]) }.join(window_lengths).map { tuple(it[1], it[2], it[3]) }
    // tree input now a tuple with 3 entries (alignment_file, model_str, relative_length)

    // choose tree reconstruction program that should be used
    def tree
    if (params.phylo_method == "iqtree2") {
        // IQTREE2 module call - will update later if module needs unique_outdir param
        tree = IQTREE2(tree_input, num_windows)
        tree.iqtree.collectFile(storeDir: "${unique_outdir}/tree_reconstruction_logs/iqtree_files")
    }
    else if (params.phylo_method == "raxml-ng") {
        // RAXMLNG module call - will update later if module needs unique_outdir param
        tree = RAXMLNG(tree_input, num_windows)
    }
    else {
        throw new Error("Phylogenetic Program'${params.phylo_method}' not supported")
    }

    tree.log.collectFile(storeDir: "${unique_outdir}/tree_reconstruction_logs/logs")

    // ALWAYS collect unrooted trees (regardless of rooting settings)
    def unrooted_tree_collection = CollectTrees(
        tree.treefile.flatten().collect().map { tuple("unrooted_trees", it) },
        unique_outdir,
    )

    // MAD rooting for best trees (if enabled)
    if (params.mad_rooting) {
        // Create proper channels for alignment and tree pairing
        def alignment_channel = alignment_windows
            .flatten()
            .map { file ->
                def baseName = file.baseName.replaceAll('\\.fasta$', '')
                tuple(baseName, file)
            }

        def tree_channel = tree.treefile
            .flatten()
            .map { file ->
                def baseName = file.baseName.replaceAll('\\.treefile$', '').replaceAll('\\.fasta$', '')
                tuple(baseName, file)
            }

        // Join alignments and trees by their base name, then create proper input tuple
        def mad_input = alignment_channel
            .join(tree_channel)
            .map { window_name, alignment_file, tree_file ->
                tuple(window_name, tree_file, alignment_file)
            }

        rooted_trees_mad = MAD_ROOTING(mad_input, unique_outdir)
        rooted_tree_collection = COLLECT_ROOTED_TREES(
            rooted_trees_mad.rooted_trees.map { it[1] }.collect(),
            rooted_trees_mad.logs.collect(),
            params.output_format,
            unique_outdir,
            main_results_dir,
        )
    }
    // Note: Unrooted trees are collected above regardless of rooting settings

    // Backbone midpoint rooting (if enabled)
    if (params.backbone_midpoint_rooting) {
        backbone_tree = GENERATE_BACKBONE_TREE(inputFile, unique_outdir)

        // Create proper channels for alignment and tree pairing (similar to MAD rooting)
        def alignment_channel = alignment_windows
            .flatten()
            .map { file ->
                def baseName = file.baseName.replaceAll('\\.fasta$', '')
                tuple(baseName, file)
            }

        def tree_channel = tree.treefile
            .flatten()
            .map { file ->
                def baseName = file.baseName.replaceAll('\\.treefile$', '').replaceAll('\\.fasta$', '')
                tuple(baseName, file)
            }

        // Join alignments and trees by their base name, then create proper input tuple
        def backbone_input = alignment_channel
            .join(tree_channel)
            .map { window_name, alignment_file, tree_file ->
                tuple(window_name, tree_file, alignment_file)
            }

        backbone_rooted_trees = BACKBONE_MIDPOINT_ROOTING(
            backbone_input,
            backbone_tree.backbone_tree,
            unique_outdir,
        )

        backbone_collection = COLLECT_BACKBONE_MIDPOINT_TREES(
            backbone_rooted_trees.rooted_trees.map { it[1] }.collect(),
            backbone_rooted_trees.constraint_trees.map { it[1] }.collect(),
            backbone_rooted_trees.logs.collect(),
            params.output_format,
            unique_outdir,
            main_results_dir,
        )
    }
}
