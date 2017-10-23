#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "cle_annotated_vcf_filter"
arguments: [
    "/usr/bin/perl", "/usr/bin/docm_and_coding_indel_selection.pl",
    $(inputs.annotated_vcf.path), $(runtime.outdir)
]
inputs:
    annotated_vcf:
        type: File
    filter:
        type: boolean
        inputBinding:
            prefix: "filter"
            position: 1
outputs:
    annotated_filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated_filtered.vcf"
