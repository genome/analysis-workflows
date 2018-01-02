#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Kallisto: Quant"
baseCommand: ["/usr/bin/kallisto"] 
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 8
arguments: [
    "quant",
    "-t", $(runtime.cores),
    "-b", "100",
    "-o", "kallisto"
    "--fr-stranded",
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
