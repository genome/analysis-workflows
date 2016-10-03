#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "vcf merge"
baseCommand: ["/usr/bin/bcftools", "concat"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.3.1"
arguments:
    - "--allow-overlaps"
    - "--remove-duplicates"
    - "--output-type"
    - "z"
    - "-o"
    - { valueFrom: $(runtime.outdir)/merged.vcf }
inputs:
    vcfs:
        type: File[]
        inputBinding:
            position: 1
        secondaryFiles: .tbi
outputs:
    merged_vcf:
        type: File
        outputBinding:
            glob: "merged.vcf"
