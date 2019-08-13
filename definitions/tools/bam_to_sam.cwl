#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'BAM to SAM conversion'
baseCommand: ["/opt/samtools/bin/samtools", "view", "-h", "-o"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
    - class: ResourceRequirement
      ramMin: 16000

arguments: [out.sam]

inputs:
    bam:
        type: File
        inputBinding:
            position: 1
outputs:
    final_sam:
        type: File
        outputBinding:
            glob: "out.sam"
