subworkflow split_normalize_swf:
    workdir:
            "."
    snakefile:
            "nodes/utils/split_normalize/Snakefile"

subworkflow sequence_length_distribution_swf:
    workdir:
            "."
    snakefile:
            "nodes/plots/sequence_length_distribution/Snakefile"

subworkflow generate_multiple_sequence_alignment_swf:
    workdir:
            "."
    snakefile:
            "nodes/utils/multiple_sequence_alignment/Snakefile"

subworkflow secondary_structure_profile_swf:
    workdir:
            "."
    snakefile:
            "nodes/utils/secondary_structure_profile/Snakefile"

subworkflow disorder_swf:
    workdir:
           "."
    snakefile:
           "nodes/encodings/disorder/Snakefile"

subworkflow disorderb_swf:
    workdir:
           "."
    snakefile:
           "nodes/encodings/disorderb/Snakefile"

subworkflow disorderc_swf:
    workdir:
           "."
    snakefile:
           "nodes/encodings/disorderc/Snakefile"
    configfile:
           "nodes/encodings/aac/config.yaml"

subworkflow aac_swf:
    workdir:
           "."
    snakefile:
           "nodes/encodings/aac/Snakefile"
    configfile:
           "nodes/encodings/aac/config.yaml"

