#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'sort BAM by position'
baseCommand: ["/usr/bin/sambamba", "sort"]
arguments: 
    ["-t", { valueFrom: $(runtime.cores) },
    "-m", "12G",
    "-o", { valueFrom: $(runtime.outdir)/PositionSorted.bam }]
requirements:
    - class: ResourceRequirement
      ramMin: 12000
      coresMin: 8
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
outputs:
    position_sorted_bam:
        type: File
        outputBinding:
            glob: "PositionSorted.bam"
