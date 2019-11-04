#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Insert an arbitrary string at a specific line of a gzipped file"
baseCommand: ["zcat"]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
    - class: ResourceRequirement
      ramMin: 4000
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
stdout: "$(inputs.input_file.basename).commented.gz"
outputs:
    output_file:
        type: stdout
