#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["/opt/samtools/bin/samtools", "sort", "-n"]
requirements:
    - class: ResourceRequirement
      ramMin: 26000
      coresMin: 8
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"

arguments:
  - prefix: -o
    valueFrom: $(inputs.cram.nameroot).NameSorted.cram
  - prefix: -@
    valueFrom: $(runtime.cores)

inputs:
  reference:
    type:
        - string
        - File
    secondaryFiles: [.fai]
    inputBinding:
        position: 1
        prefix: '--reference'
  cram:
    type: File
    inputBinding:
        position: 2

outputs:
  name_sorted_cram:
    type: File
    outputBinding:
      glob: $(inputs.cram.nameroot).NameSorted.cram

