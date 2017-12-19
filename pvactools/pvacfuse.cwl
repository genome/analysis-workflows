#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run pVACfuse"

baseCommand: [
    "pvacfuse",
    "run"
]
requirements:
    - class: DockerRequirement
      dockerPull: registry.gsc.wustl.edu/ssiebert/docker-pvactools
arguments:
    - position: 5
      valueFrom: $(runtime.outdir)
inputs:
    input_file:
        type: File
        inputBinding:
            position: 1
    sample_name:
        type: string
        inputBinding:
            position: 2
    alleles:
        type: string[]
        inputBinding:
            position: 3
            itemSeparator: ','
            separate: false
            prefix: ""
    prediction_algorithms:
        type: string[]
        inputBinding:
            position: 4
    epitope_lengths:
        type: int[]?
        inputBinding:
            prefix: "-e"
            itemSeparator: ','
            separate: false
    iedb_install_directory:
        type: string?
        inputBinding:
            prefix: "--iedb-install-directory"
    binding_threshold:
        type: int?
        inputBinding:
            prefix: "-b"
    iedb_retries:
        type: int?
        inputBinding:
            prefix: "-r"
    keep_tmp_files:
        type: boolean?
        inputBinding:
            prefix: "-k"
    peptide_sequence_length:
        type: int?
        inputBinding:
            prefix: "-l"
        default: 21
    net_chop_method:
        type:
            - "null"
            - type: enum
              symbols: ["cterm", "20s"]
        inputBinding:
            prefix: "--net-chop-method"
    netmhc_stab:
        type: boolean?
        inputBinding:
            prefix: "--netmhc-stab"
    top_result_per_mutation:
        type: boolean?
        inputBinding:
            prefix: "-t"
    top_score_metric:
        type:
            - "null"
            - type: enum
              symbols: ["lowest", "median"]
        inputBinding:
            prefix: "-m"
    net_chop_threshold:
        type: float?
        inputBinding:
            prefix: "--net-chop-threshold"
    additional_report_columns:
        type:
            - "null"
            - type: enum
              symbols: ["sample_name"]
        inputBinding:
            prefix: "-a"
    fasta_size:
        type: int?
        inputBinding:
            prefix: "-s"
    downstream_sequence_length:
        type: string?
        inputBinding:
            prefix: "-d"
    exclude_nas:
        type: boolean?
        inputBinding:
            prefix: "--exclude-NAs"
outputs:
    mhc_i_unfiltered_neoepitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_I/*.combined.parsed.tsv"
    mhc_i_filtered_neoepitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_I/*.final.tsv"
    mhc_ii_unfiltered_neoepitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_II/*.combined.parsed.tsv"
    mhc_ii_filtered_neoepitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_II/*.final.tsv"
