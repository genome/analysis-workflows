#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "annotated_vcf_filter"
arguments: [
    "/usr/bin/perl", "/usr/bin/docm_filter.pl",
    $(inputs.annotated_vcf.path), $(runtime.outdir)
]
inputs:
    annotated_vcf:
        type: File
outputs:
    annotated_filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated_filtered.vcf"
