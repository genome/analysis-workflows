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
      dockerPull: "mgibio/samtools-cwl:1.0.0"
arguments: ["$(inputs.name).merged.bam"]
inputs:
    bams:
        type: File[]
        inputBinding:
            position: 1
    name:
        type: string
outputs:
    merged_bam:
        type: File
        outputBinding:
            glob: "$(inputs.name).merged.bam"

