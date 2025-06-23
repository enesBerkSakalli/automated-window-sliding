//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml

class Utils {

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by required channel priority.
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that they are in the right order
        def channel_priority_violation = false
        def n = required_channels_in_order.size()
        for (int i = 0; i < n - 1; i++) {
            channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
        }

        if (channels_missing | channel_priority_violation) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

    //
    // Check resource limits against maximum allowed values
    //
    public static Object check_max(obj, type) {
        if (type == 'memory') {
            try {
                if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                    return params.max_memory as nextflow.util.MemoryUnit
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'time') {
            try {
                if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                    return params.max_time as nextflow.util.Duration
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'cpus') {
            try {
                return Math.min( obj, params.max_cpus as int )
            } catch (all) {
                println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
                return obj
            }
        }
    }
}
