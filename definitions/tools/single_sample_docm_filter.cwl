#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "docm filter"
arguments: [
    "/usr/bin/perl", "/usr/bin/single_sample_docm_filter.pl",
    $(inputs.docm_out.path), $(runtime.outdir)
]
inputs:
    docm_out:
        type: File
outputs:
    docm_filter_out:
        type: File
        outputBinding:
            glob: "docm_filter_out.vcf"
