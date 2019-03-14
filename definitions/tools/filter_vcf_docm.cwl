#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Filter variants from the DoCM detector"
baseCommand: ["/usr/bin/perl", "/usr/bin/docm_filter.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
arguments: [
    $(runtime.outdir)/docm_filtered_variants.vcf
]
inputs:
    docm_raw_variants:
        type: File
        inputBinding:
            position: -4
    normal_bam:
        type: File
        inputBinding:
            position: -3
    tumor_bam:
        type: File
        inputBinding:
            position: -2
    filter_docm_variants:
        type: int
        inputBinding:
            position: -1
outputs:
    docm_filtered_variants:
        type: File
        outputBinding:
            glob: "docm_filtered_variants.vcf"
