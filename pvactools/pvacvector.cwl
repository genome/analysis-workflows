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
    input_vcf:
        type: File?
        inputBinding:
            prefix: "-v"
    input_n_mer:
        type: int?
        inputBinding:
            prefix: "-n"
        default: 21
outputs:
    vector_fasta:
        type: File
        outputBinding:
            glob: "Test_results.fa"
    vector_jpg:
        type: File?
        outputBinding:
            glob: "vector.jpg"
