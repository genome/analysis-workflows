#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Samtools: merge"
baseCommand: ["/opt/samtools/bin/samtools", "merge"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
arguments: ["merged.bam"]
inputs:
    bams:
        type: File[]
        inputBinding:
            position: 1
outputs:
    merged_bam:
        type: File
        outputBinding:
            glob: "merged.bam"

