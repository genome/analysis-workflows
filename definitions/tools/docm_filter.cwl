#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "docm filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/docm_filter.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
arguments: [
    $(inputs.docm_out.path), $(inputs.normal_cram.path), $(inputs.tumor_cram.path), $(runtime.outdir)
]
inputs:
    docm_out:
        type: File
    normal_cram:
        type: File
    tumor_cram:
        type: File
outputs:
    docm_filter_out:
        type: File
        outputBinding:
            glob: "docm_filter_out.vcf"
