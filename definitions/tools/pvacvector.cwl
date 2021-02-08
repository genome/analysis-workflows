#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run pVACvector"

baseCommand: [
    "ln", "-s"
]
arguments: [
    { valueFrom: "$TMPDIR", shellQuote: false },
    "/tmp/pvacseq",
    { valueFrom: " && ", shellQuote: false },
    "export", "TMPDIR=/tmp/pvacseq",
    { valueFrom: " && ", shellQuote: false },
    "/usr/local/bin/pvacvector",
    "run",
    "--iedb-install-directory", "/opt/iedb",
    { position: 5, valueFrom: $(runtime.outdir) },
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "griffithlab/pvactools:2.0.0"
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: $(inputs.n_threads)
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
    top_score_metric:
        type:
            - "null"
            - type: enum
              symbols: ["lowest", "median"]
        inputBinding:
            prefix: "-m"
    iedb_retries:
        type: int?
        inputBinding:
            prefix: "-r"
    keep_tmp_files:
        type: boolean?
        inputBinding:
            prefix: "-k"
    input_vcf:
        type: File?
        inputBinding:
            prefix: "-v"
    input_n_mer:
        type: int?
        inputBinding:
            prefix: "-n"
        default: 21
    n_threads:
        type: int?
        inputBinding:
            prefix: "--n-threads"
        default: 8
    spacers:
        type: string?
        inputBinding:
            prefix: "--spacers"
        default: "None,AAY,HHHH,GGS,GPGPG,HHAA,AAL,HH,HHC,HHH,HHHD,HHL,HHHC"
    max_clip_length:
        type: int?
        inputBinding:
            prefix: "--max-clip-length"
        default: 3
outputs:
    vector_fasta:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name)_results.fa"
    vector_dna_fasta:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name)_results.dna.fa"
    vector_jpg:
        type: File?
        outputBinding:
            glob: "vector.jpg"
