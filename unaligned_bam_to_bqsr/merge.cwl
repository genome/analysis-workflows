#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'merge BAMs'
baseCommand: ["/usr/local/bin/samtools", "merge"]
arguments: ["AlignedMerged.bam"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/samtools-1.3.1-2:2"
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
