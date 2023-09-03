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
    { position: 5, valueFrom: "pvacfuse_predictions" },
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "griffithlab/pvactools:4.0.0"
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: $(inputs.n_threads)
inputs:
    input_fusions:
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
    allele_specific_binding_threshold:
        type: boolean?
        inputBinding:
            prefix: "--allele-specific-binding-thresholds"
    aggregate_inclusion_binding_threshold:
        type: int?
        inputBinding:
            prefix: "--aggregate-inclusion-binding-threshold"
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
    problematic_amino_acids:
        type: string[]?
        inputBinding:
            prefix: "--problematic-amino-acids"
            itemSeparator: ','
    net_chop_threshold:
        type: float?
        inputBinding:
            prefix: "--net-chop-threshold"
    run_reference_proteome_similarity:
        type: boolean?
        inputBinding:
            prefix: "--run-reference-proteome-similarity"
    peptide_fasta:
        type: File?
        inputBinding:
            prefix: "--peptide-fasta"
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
    star_fusion_file:
        type: File?
        inputBinding:
            prefix: "--starfusion-file"
    read_support:
        type: int?
        inputBinding:
            prefix: "--read-support"
    expn_val:
        type: float?
        inputBinding:
            prefix: "--expn-val"
outputs:
    mhc_i_all_epitopes:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/MHC_Class_I/$(inputs.sample_name).all_epitopes.tsv"
    mhc_i_aggregated_report:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/MHC_Class_I/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    mhc_i_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/MHC_Class_I/$(inputs.sample_name).filtered.tsv"
    mhc_ii_all_epitopes:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/MHC_Class_II/$(inputs.sample_name).all_epitopes.tsv"
    mhc_ii_aggregated_report:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/MHC_Class_II/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    mhc_ii_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/MHC_Class_II/$(inputs.sample_name).filtered.tsv"
    combined_all_epitopes:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/combined/$(inputs.sample_name).all_epitopes.tsv"
    combined_aggregated_report:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/combined/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    combined_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "pvacfuse_predictions/combined/$(inputs.sample_name).filtered.tsv"
    pvacfuse_predictions:
        type: Directory
        outputBinding:
            glob: "pvacfuse_predictions"
