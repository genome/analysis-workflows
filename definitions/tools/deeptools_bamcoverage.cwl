#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "deepTools bamCoverage"
baseCommand: ['/usr/local/bin/bamCoverage']
requirements:
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/deeptools:3.1.3--py35h470a237_0"
arguments:
    ['-o', "$(inputs.bam.basename).$(inputs.outfile_format)", '-p', $(runtime.cores)]
inputs:
    bam:
        type: File
        inputBinding:
            prefix: '-b'
    outfile_format:
        type:
            type: enum
            symbols: ["bigwig", "bedgraph"]
        default: "bigwig"
        inputBinding:
            prefix: '-of'
    scale_factor:
        type: float?
        inputBinding:
            prefix: '--scaleFactor'
    filter_rna_strand:
        type:
            - "null"
            - type: enum
              symbols: ["forward", "reverse"]
        inputBinding:
            prefix: '--filterRNAstrand'
    bin_size:
        type: int?
        inputBinding:
            prefix: '-bs'
    region:
        type: string?
        inputBinding:
            prefix: '-r'
    blacklist_file:
        type: File?
        inputBinding:
            prefix: '-bl'
    verbose:
        type: boolean?
        inputBinding:
            prefix: '-v'
    effective_genome_size:
        type: int?
        inputBinding:
            prefix: '--effectiveGenomeSize'
    normalize_using:
        type:
            type: enum
            symbols: ["RPKM", "CPM", "BPM", "RPGC", "None"]
        default: "None"
        inputBinding:
            prefix: '--normalizeUsing'
    exact_scaling:
        type: boolean?
        inputBinding:
            prefix: '--exactScaling'
    ignore_for_normalization:
        type:
            - "null"
            - type: array
              items: string
        inputBinding:
            prefix: '-ignore'
            itemSeparator: ' '
    skip_non_covered_regions:
        type: boolean?
        inputBinding:
            prefix: '--skipNAs'
    smooth_length:
        type: int?
        inputBinding:
            prefix: '--smoothLength'
    extend_reads:
        type: boolean?
        inputBinding:
            prefix: '-e'
    ignore_duplicates:
        type: boolean?
        inputBinding:
            prefix: '--ignoreDuplicates'
    min_mapping_quality:
        type: int?
        inputBinding:
            prefix: '--minMappingQuality'
    center_reads:
        type: boolean?
        inputBinding:
            prefix: '--centerReads'
    sam_flag_include:
        type: int?
        inputBinding:
            prefix: '--samFlagInclude'
    sam_flag_exclude:
        type: int?
        inputBinding:
            prefix: '--samFlagExclude'
    min_fragment_length:
        type: int?
        inputBinding:
            prefix: '--minFragmentLength'
    max_fragment_length:
        type: int?
        inputBinding:
            prefix: '--maxFragmentLength'
outputs:
    outfile:
        type: File
        outputBinding:
            glob: "$(inputs.bam.basename).$(inputs.outfile_format)"
