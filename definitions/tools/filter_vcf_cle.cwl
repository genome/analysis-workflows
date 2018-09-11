#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "cle_annotated_vcf_filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/docm_and_coding_indel_selection.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
arguments: [$(inputs.vcf.path), $(runtime.outdir)]
inputs:
    vcf:
        type: File
    filter:
        type: boolean
        inputBinding:
            prefix: "filter"
            position: 1
outputs:
    cle_filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated_filtered.vcf"
