#!/usr/bin/env cwl-runner

class: CommandLineTool

cwlVersion: v1.0

baseCommand: [ "samtools", "view" ]

requirements:
  - class: DockerRequirement
    dockerPull: zlskidmore/samtools:1.9
  - class: ResourceRequirement
    ramMin: 20000
    coresMin: 8

arguments:
  - valueFrom: "8"
    position: 1
    prefix: "--threads"
  - valueFrom: "file.bam"
    position: 2
    prefix: "-o"
  - valueFrom: "-b"
    position: 3

inputs:
  file:
    type: File
    inputBinding:
      position: 4

outputs:
  bam_file:
    type: File
    outputBinding:
      glob: file.bam
