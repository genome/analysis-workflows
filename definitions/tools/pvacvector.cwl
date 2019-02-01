#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run pVACvector"

baseCommand: [
    "pvacvector",
    "run"
]
requirements:
    - class: DockerRequirement
      dockerPull: "griffithlab/pvactools:1.3.0"
arguments:
    - position: 5
      valueFrom: $(runtime.outdir)
    - position: 6
      valueFrom: "--iedb-install-directory"
    - position: 7
      valueFrom: "/opt/iedb"
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
    binding_threshold:
        type: int?
        inputBinding:
            prefix: "-b"
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
outputs:
    vector_fasta:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name)_results.fa"
    vector_jpg:
        type: File?
        outputBinding:
            glob: "vector.jpg"
