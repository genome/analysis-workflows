#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools sort"
baseCommand: ["/opt/samtools/bin/samtools", "sort"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
      coresMin: 1
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"

arguments:
  - prefix: -o
    valueFrom: $(runtime.outdir)/$(inputs.output_filename)
  - prefix: -@
    valueFrom: $(runtime.cores)

inputs:
  output_filename:
    type: string
    default: sorted.bam
  input_bam:
    type: File
    inputBinding:
      position: 50

outputs:
  sorted_bam:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

