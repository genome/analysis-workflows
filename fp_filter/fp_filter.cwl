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
    "--bam-readcount", "/usr/bin/bam-readcount",
    "--samtools", "/opt/samtools/bin/samtools",
    "--output", { valueFrom: $(runtime.outdir)/$(inputs.output_vcf_basename).vcf }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "--reference"
            position: 1
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
    output_vcf_basename:
        type: string?
        default: fpfilter
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_basename).vcf

