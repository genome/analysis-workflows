#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Sambamba: merge"
baseCommand: ["/usr/bin/sambamba_merge"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
      coresMin: 4
arguments: [
    "$(runtime.cores)",
    "$(runtime.outdir)/merged.bam"
]
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
