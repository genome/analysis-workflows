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
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
arguments: ["$(inputs.name).merged.cram", { prefix: "--threads", valueFrom: $(runtime.cores) }, { prefix: '-O', valueFrom: "CRAM"}]
inputs:
    crams:
        type: File[]
        inputBinding:
            position: 1
    name:
        type: string
    reference:
        type:
            - string
            - File
        inputBinding:
            position: 2
            prefix: '--reference'

outputs:
    merged_cram:
        type: File
        outputBinding:
            glob: "$(inputs.name).merged.cram"

