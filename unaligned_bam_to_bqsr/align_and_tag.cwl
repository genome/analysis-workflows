#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "align with bwa_mem and tag"

baseCommand: ["/bin/bash", "/usr/bin/alignment_helper.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 16000
stdout: "refAlign.bam"
arguments:
    - position: 5
      valueFrom: $(runtime.cores)
inputs:
    readgroup:
        type: string
        inputBinding:
            position: 1
    reference:
        type: string
        inputBinding:
            position: 2
    fastq:
        type: File
        inputBinding:
            position: 3
    fastq2:
        type: File
        inputBinding:
            position: 4
outputs:
    aligned_bam:
            type: stdout
