#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools index"
baseCommand: ["index"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/samtools:1.3.1"
arguments:
    - position: 2
      valueFrom: $(runtime.outdir)/$(inputs.bam.basename).bai
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
outputs:
    bam_index:
        type: File
        outputBinding:
            glob: $(inputs.bam.basename).bai
