#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'BAM to CRAM conversion'
baseCommand: ["/opt/samtools/bin/samtools", "view", "-C"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
stdout: "final.cram"
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
outputs:
    cram:
        type: stdout
