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
    "/usr/local/bin/pvacseq", "run",
    "--iedb-install-directory", "/opt/iedb",
    "--blastp-path", "/opt/ncbi-blast-2.12.0+/bin/blastp",
    "--pass-only",
    { position: 5, valueFrom: "pvacseq_predictions" },
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "griffithlab/pvactools:3.1.1"
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: $(inputs.n_threads)
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
    run_reference_proteome_similarity:
        type: boolean?
        inputBinding:
            prefix: "--run-reference-proteome-similarity"
    blastp_db:
        type:
            - "null"
            - type: enum
              symbols: ["refseq_select_prot", "refseq_protein"]
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
    tumor_purity:
        type: float?
        inputBinding:
            prefix: "--tumor-purity"
    n_threads:
        type: int?
        inputBinding:
            prefix: "--n-threads"
        default: 8
outputs:
    mhc_i_all_epitopes:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/MHC_Class_I/$(inputs.sample_name).all_epitopes.tsv"
    mhc_i_aggregated_report:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/MHC_Class_I/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    mhc_i_aggregated_metrics_file:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/MHC_Class_I/$(inputs.sample_name).all_epitopes.aggregated.metrics.json"
    mhc_i_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/MHC_Class_I/$(inputs.sample_name).filtered.tsv"
    mhc_ii_all_epitopes:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/MHC_Class_II/$(inputs.sample_name).all_epitopes.tsv"
    mhc_ii_aggregated_report:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/MHC_Class_II/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    mhc_ii_aggregated_metrics_file:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/MHC_Class_II/$(inputs.sample_name).all_epitopes.aggregated.metrics.json"
    mhc_ii_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/MHC_Class_II/$(inputs.sample_name).filtered.tsv"
    combined_all_epitopes:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/combined/$(inputs.sample_name).all_epitopes.tsv"
    combined_aggregated_report:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/combined/$(inputs.sample_name).all_epitopes.aggregated.tsv"
    combined_aggregated_metrics_file:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/combined/$(inputs.sample_name).all_epitopes.aggregated.metrics.json"
    combined_filtered_epitopes:
        type: File?
        outputBinding:
            glob: "pvacseq_predictions/combined/$(inputs.sample_name).filtered.tsv"
    pvacseq_predictions:
        type: Directory
        outputBinding:
            glob: "pvacseq_predictions"
