#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "merge VCF files from non-overlapping sample sets"
baseCommand: ["/opt/bcftools/bin/bcftools", "merge"]
requirements:
    - class: DockerRequirement
      dockerPull: mgibio/bcftools-cwl:1.3.1
    - class: ResourceRequirement
      ramMin: 4000
arguments:
    - "-o"
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
