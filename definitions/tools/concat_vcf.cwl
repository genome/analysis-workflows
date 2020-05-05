#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "concatenate or combine multiple VCF files that contain variants from the same samples"
baseCommand: ["/opt/bcftools/bin/bcftools", "concat"]
requirements:
    - class: DockerRequirement
      dockerPull: mgibio/bcftools-cwl:1.3.1
    - class: ResourceRequirement
      ramMin: 4000
arguments:
    - "--allow-overlaps"
    - "--remove-duplicates"
    - "--output-type"
    - "z"
    - "-o"
    - { valueFrom: $(runtime.outdir)/$(inputs.merged_vcf_basename).vcf.gz }
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
            glob: $(inputs.merged_vcf_basename).vcf.gz