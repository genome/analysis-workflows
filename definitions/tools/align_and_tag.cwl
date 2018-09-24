#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "align with bwa_mem and tag"

baseCommand: ["/bin/bash", "/usr/bin/alignment_helper.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 20000
    - class: DockerRequirement
      dockerPull: "mgibio/alignment_helper-cwl:1.0.0"
stdout: "refAlign.bam"
arguments:
    - position: 4
      valueFrom: $(runtime.cores)
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
    readgroup:
        type: string
        inputBinding:
            position: 2
    reference:
        type: string
        inputBinding:
            position: 3
outputs:
    aligned_bam:
            type: stdout
