#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Insert an arbitrary string at a specific line of a gzipped file"
baseCommand: ["zcat"]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    - class: ResourceRequirement
      ramMin: 4000
    - class: StepInputExpressionRequirement
arguments:
    [ { valueFrom: $(inputs.input_file.path) } , { shellQuote: false, valueFrom: "|" },
      "awk", "-v", { valueFrom: n=$(inputs.line_number) }, "-v", { valueFrom: s=$(inputs.some_text) }, 'NR == n {print s} {print}', { shellQuote: false, valueFrom: "|" },
      "/opt/htslib/bin/bgzip"
    ]
inputs:
    input_file:
        type: File
    line_number:
        type: int
    some_text:
        type: string
    output_name:
        type: string?
        default: "$(inputs.input_file.basename).commented.gz"
stdout: $(inputs.output_name)
outputs:
    output_file:
        type: stdout
