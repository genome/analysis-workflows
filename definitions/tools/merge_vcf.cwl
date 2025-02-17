#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "vcf merge"
baseCommand: ["/opt/bcftools/bin/bcftools", "merge"]
requirements:
    - class: DockerRequirement
      dockerPull: mgibio/bcftools-cwl:1.12
    - class: ResourceRequirement
      ramMin: 4000
arguments:
    - "--filter-logic"
    - "x"
    - "--merge"
    - "none"
    - "--output"
    - { valueFrom: $(runtime.outdir)/$(inputs.merged_vcf_basename).vcf }
inputs:
    vcfs:
        type: File[]
        inputBinding:
            position: 1
        secondaryFiles: [.tbi]
    merged_vcf_basename:
        type: string?
        default: 'merged'
outputs:
    merged_vcf:
        type: File
        outputBinding:
            glob: $(inputs.merged_vcf_basename).vcf
