#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "docm filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/docm_filter.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.3.1"
    - class: ResourceRequirement
      ramMin: 4000
arguments: [
    $(inputs.docm_out.path), $(inputs.normal_bam.path), $(inputs.tumor_bam.path), $(runtime.outdir)
]
inputs:
    docm_out:
        type: File
    normal_bam:
        type: File
    tumor_bam:
        type: File
outputs:
    docm_filter_out:
        type: File
        outputBinding:
            glob: "docm_filter_out.vcf"
