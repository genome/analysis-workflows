#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run pVACfuse"

baseCommand: [
    "ln", "-s"
]
arguments: [
    { valueFrom: "$TMPDIR", shellQuote: false },
    "/tmp/pvacseq",
    { valueFrom: " && ", shellQuote: false },
    "export", "TMPDIR=/tmp/pvacseq",
    { valueFrom: " && ", shellQuote: false },
    "/usr/local/bin/pvacfuse",
    "run",
    "--iedb-install-directory", "/opt/iedb",
    { position: 5, valueFrom: $(runtime.outdir) },
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "griffithlab/pvactools:2.0.1"
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: $(inputs.n_threads)
inputs:
    input_file:
        type:
            - File
            - Directory
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
    epitope_lengths_class_i:
        type: int[]?
        inputBinding:
            prefix: "-e1"
            itemSeparator: ','
    epitope_lengths_class_ii:
        type: int[]?
        inputBinding:
            prefix: "-e2"
            itemSeparator: ','
    binding_threshold:
        type: int?
        inputBinding:
            prefix: "-b"
    percentile_threshold:
        type: int?
        inputBinding:
            prefix: "--percentile-threshold"
    iedb_retries:
        type: int?
        inputBinding:
            prefix: "-r"
    keep_tmp_files:
        type: boolean?
        inputBinding:
            prefix: "-k"
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
    run_reference_proteome_similarity:
        type: boolean?
        inputBinding:
            prefix: "--run-reference-proteome-similarity"
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
    n_threads:
        type: int?
        inputBinding:
            prefix: "--n-threads"
        default: 8
outputs:
    mhc_i_all_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_I/$(inputs.sample_name).all_epitopes.tsv"
    mhc_i_aggregated_report:
        type: File?
        outputBinding:
            glob: "MHC_Class_I/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    mhc_i_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_I/$(inputs.sample_name).filtered.tsv"
    mhc_ii_all_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_II/$(inputs.sample_name).all_epitopes.tsv"
    mhc_ii_aggregated_report:
        type: File?
        outputBinding:
            glob: "MHC_Class_II/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    mhc_ii_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_II/$(inputs.sample_name).filtered.tsv"
    combined_all_epitopes:
        type: File?
        outputBinding:
            glob: "combined/$(inputs.sample_name).all_epitopes.tsv"
    combined_aggregated_report:
        type: File?
        outputBinding:
            glob: "combined/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    combined_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "combined/$(inputs.sample_name).filtered.tsv"
