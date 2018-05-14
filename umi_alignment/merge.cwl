#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'merge BAMs'
baseCommand: ["/opt/samtools/bin/samtools", "merge"]
arguments: ["AlignedMerged.bam"]
inputs:
    bams:
        type: File[]
        inputBinding:
            position: 1
outputs:
    merged_bam:
        type: File
        outputBinding:
            glob: "AlignedMerged.bam"
