#!/usr/bin/env cwl-runner

class: CommandLineTool

cwlVersion: v1.0

baseCommand: [ "samtools", "index" ]

requirements:
  - class: DockerRequirement
    dockerImageId: zlskidmore/samtools:1.9
    dockerPull: zlskidmore/samtools:1.9
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.bam_file)
  - class: ResourceRequirement
    ramMin: 4000

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1

outputs:
  input_index:
    type: File
    secondaryFiles: [.bai, .crai]
    outputBinding:
      glob: "*.bam"
