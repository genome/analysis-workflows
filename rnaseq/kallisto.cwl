#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Kallisto: Quant"
baseCommand: ["/opt/kallisto/kallisto_linux-v0.43.1/kallisto"] #TODO - switch back to /usr/bin/kallisto after rebuilding container
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 8
arguments: [
    "quant",
    "-t", $(runtime.cores),
    "-b", "100",
    "-o", "kallisto",
    $(inputs.fastq1),
    $(inputs.fastq2)
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
            glob: "kallisto/abundance.tsv"
