#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Kallisto: Quant"
baseCommand: ["/opt/kallisto/kallisto_linux-v0.43.1/kallisto"] #TODO - switch back to /usr/bin/kallisto after rebuilding container
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
    - class: MultipleInputFeatureRequirement
      ramMin: 32000
      coresMin: 8
arguments: [
    "quant",
    "-t", $(runtime.cores),
    "-b", "100",
    "-o", "kallisto"
]
inputs:
    kallisto_index:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    fastqs:
        type:
            type: array
            items: 
                type: array
                items: File
        inputBinding:
            position: 2
outputs:
    expression_transcript_table:
        type: File
        outputBinding:
            glob: "kallisto/abundance.tsv"
    expression_transcript_h5:
        type: File
        outputBinding:
            glob: "kallisto/abundance.h5"
