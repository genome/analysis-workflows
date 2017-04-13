#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "False Positive filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/fpfilter.pl"]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
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

