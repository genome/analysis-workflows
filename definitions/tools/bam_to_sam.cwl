#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'BAM to SAM conversion'
baseCommand: ["/usr/local/bin/samtools", "view", "-h", "-o"]
requirements:
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    - class: ResourceRequirement
      ramMin: 16000

arguments: ["$(inputs.bam.nameroot).sam"]

inputs:
    bam:
        type: File
        inputBinding:
            position: 1
outputs:
    final_sam:
        type: File
        outputBinding:
            glob: "$(inputs.bam.nameroot).sam"
