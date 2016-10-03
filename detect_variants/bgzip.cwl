#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "bgzip VCF"
baseCommand: ["/opt/samtools/bin/bgzip"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/tabix-cwl:1.3.1"
stdout: $(inputs.file.basename).gz
arguments:
    ["-c"]
inputs:
    file:
        type: File
        inputBinding:
            position: 1
outputs:
    bgzipped_file:
        type: stdout

