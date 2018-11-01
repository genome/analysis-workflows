#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'BAM to CRAM conversion'
baseCommand: ["/opt/samtools/bin/samtools", "view", "-C"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
stdout: $(inputs.name)
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-T"
            position: 1
    bam:
        type: File
        inputBinding:
            position: 2
    name:
        type: string?
        default: 'final.cram'
outputs:
    cram:
        type: stdout
