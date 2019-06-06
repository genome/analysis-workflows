#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Insert an arbitrary string at a specific line of a file"
baseCommand: ["awk"]
requirements:
    - class: DockerRequirement
      dockerPull: 'ubuntu:xenial'
arguments:
    [ "-v", { valueFrom: n=$(inputs.line_number) }, "-v", { valueFrom: s=$(inputs.some_text) }, 'NR == n {print s} {print}']
inputs:
    input_file:
        type: File
        inputBinding:
            position: 7
    line_number:
        type: int
    some_text:
        type: string
stdout: "$(inputs.input_file.basename).commented"
outputs:
    output_file:
        type: stdout
