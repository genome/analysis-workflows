#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "vcf merge"
baseCommand: ["/opt/bcftools/bin/bcftools", "concat"]
arguments:
    - "--allow-overlaps"
    - "--remove-duplicates"
    - "--output-type"
    - "z"
    - "-o"
    - { valueFrom: $(runtime.outdir)/merged.vcf.gz }
inputs:
    vcfs:
        type: File[]
        inputBinding:
            position: 1
        secondaryFiles: [.tbi]
outputs:
    merged_vcf:
        type: File
        outputBinding:
            glob: "merged.vcf.gz"
