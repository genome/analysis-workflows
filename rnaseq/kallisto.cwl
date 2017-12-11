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
    "-t", $(runtime.cores),
    "-b", "100",
    "-o", "transcriptQuant.tsv"
    "$(inputs.fastqs)
]
inputs:
    kallisto_index:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    fastq1:
        type: File
        inputBinding:
        position: 2
    fastq2:
        type: File
        inputBinding:        
        position: 3

outputs:
    expression_transcript_table:
        type: File
        outputBinding:
            glob: "transcriptQuant.tsv"
