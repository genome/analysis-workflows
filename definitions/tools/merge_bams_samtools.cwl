#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Samtools: merge"
baseCommand: ["/usr/local/bin/samtools", "merge"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/samtools:0.1.19--hfb9b9cc_8"
arguments: ["$(inputs.name).merged.bam", { prefix: "-@", valueFrom: $(runtime.cores) }]
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

