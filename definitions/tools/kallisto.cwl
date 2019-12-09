#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Kallisto: Quant"
baseCommand: ["/usr/bin/kallisto"]
requirements:
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 8
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq:1.0.0"
arguments: [
    "quant",
    "-t", $(runtime.cores),
    "-b", "100",
    "--fusion",
    "-o", "kallisto",
]
inputs:
    kallisto_index:
        type: File
        inputBinding:
            prefix: "-i"
            position: 2
    fastqs:
        type:
            type: array
            items:
                type: array
                items: File
        inputBinding:
            position: 3
    strand:
        type:
          - "null"
          - type: enum
            symbols: ["first", "second", "unstranded"]
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.strand) {
                        if (inputs.strand == 'first') {
                            return ['--fr-stranded'];
                        } else if (inputs.strand == 'second') {
                            return ['--rf-stranded'];
                        } else {
                            return [];
                        }
                    } else {
                            return [];
                    }
                }
            position: 1
outputs:
    expression_transcript_table:
        type: File
        outputBinding:
            glob: "kallisto/abundance.tsv"
    expression_transcript_h5:
        type: File
        outputBinding:
            glob: "kallisto/abundance.h5"
    fusion_evidence:
        type: File
        outputBinding:
            glob: "kallisto/fusion.txt"
