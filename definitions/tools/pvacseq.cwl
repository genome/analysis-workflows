#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run pVACseq"

baseCommand: [
    "ln", "-s"
]
arguments: [
    { valueFrom: "$TMPDIR", shellQuote: false },
    "/tmp/pvacseq",
    { valueFrom: " && ", shellQuote: false },
    "export", "TMPDIR=/tmp/pvacseq",
    { valueFrom: " && ", shellQuote: false },
    "/opt/conda/bin/pvacseq", "run",
    "--iedb-install-directory", "/opt/iedb",
    "--pass-only",
    { position: 5, valueFrom: $(runtime.outdir) },
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "griffithlab/pvactools:1.4.2"
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 8
inputs:
    input_vcf:
        type: File
        secondaryFiles: ['.tbi']
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
    binding_threshold:
        type: int?
        inputBinding:
            prefix: "-b"
    allele_specific_binding_thresholds:
        type: boolean?
        inputBinding:
            prefix: "--allele-specific-binding-thresholds"
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
    normal_sample_name:
        type: string?
        inputBinding:
            prefix: "--normal-sample-name"
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
    phased_proximal_variants_vcf:
        type: File?
        secondaryFiles: [".tbi"]
        inputBinding:
            prefix: "-p"
    minimum_fold_change:
        type: float?
        inputBinding:
            prefix: "-c"
    normal_cov:
        type: int?
        inputBinding:
            prefix: "--normal-cov"
    tdna_cov:
        type: int?
        inputBinding:
            prefix: "--tdna-cov"
    trna_cov:
        type: int?
        inputBinding:
            prefix: "--trna-cov"
    normal_vaf:
        type: float?
        inputBinding:
            prefix: "--normal-vaf"
    tdna_vaf:
        type: float?
        inputBinding:
            prefix: "--tdna-vaf"
    trna_vaf:
        type: float?
        inputBinding:
            prefix: "--trna-vaf"
    expn_val:
        type: float?
        inputBinding:
            prefix: "--expn-val"
    maximum_transcript_support_level:
        type:
            - "null"
            - type: enum
              symbols: ["1", "2", "3", "4", "5"]
        inputBinding:
            prefix: "--maximum-transcript-support-level"
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
    mhc_i_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_I/$(inputs.sample_name).filtered.tsv"
    mhc_i_ranked_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_I/$(inputs.sample_name).filtered.condensed.ranked.tsv"
    mhc_ii_all_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_II/$(inputs.sample_name).all_epitopes.tsv"
    mhc_ii_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_II/$(inputs.sample_name).filtered.tsv"
    mhc_ii_ranked_epitopes:
        type: File?
        outputBinding:
            glob: "MHC_Class_II/$(inputs.sample_name).filtered.condensed.ranked.tsv"
    combined_all_epitopes:
        type: File?
        outputBinding:
            glob: "combined/$(inputs.sample_name).all_epitopes.tsv"
    combined_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "combined/$(inputs.sample_name).filtered.tsv"
    combined_ranked_epitopes:
        type: File?
        outputBinding:
            glob: "combined/$(inputs.sample_name).filtered.condensed.ranked.tsv"
