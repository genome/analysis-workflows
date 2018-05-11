#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Bam To Cram"
baseCommand: ["/usr/bin/bam_to_cram"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 16000
      coresMin:
arguments: [
    "$(runtime.outdir)/merged.cram"
]
inputs:
    reference_index:
        type: string
        inputBinding:
            position: -2
    bam:
        type: File
        inputBinding:
            position: -1
outputs:
    cram:
        type: File
        secondaryFiles: [.crai]
        outputBinding:
            glob: "merged.cram"
