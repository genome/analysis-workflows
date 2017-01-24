#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "False Positive filter"
baseCommand: ["/usr/bin/fpfilter.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/fpfilter-cwl:v1"
arguments:
    ["--sample", "TUMOR",
    "--output", { valueFrom: $(runtime.outdir)/filtered.vcf }]
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "--reference"
            position: 1
        secondaryFiles: [.fai]
    bam:
        type: File
        inputBinding:
            prefix: "--bam-file"
            position: 2
        secondaryFiles: [.bai]
    vcf:
        type: File
        inputBinding:
            prefix: "--vcf-file"
            position: 3
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: "filtered.vcf"

