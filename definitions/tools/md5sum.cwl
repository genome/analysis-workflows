#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/md5sum']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:bionic"
    - class: ResourceRequirement
      ramMin: 16000
    - class: StepInputExpressionRequirement

inputs:
    input_files:
        type: File[]
        inputBinding:
            position: 1
    output_name:
        type: string?
        default: $(inputs.input_files[0].basename)
stdout: "$(inputs.output_name).md5"
outputs:
    md5sum:
        type: stdout
