#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add GT tags"
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.3.1"
    - class: ResourceRequirement
      ramMin: 4000
arguments: [
    "/usr/bin/perl", "/usr/bin/add_strelka_gt.pl",
    $(inputs.vcf.path), $(runtime.outdir)
]
inputs:
    vcf:
        type: File
outputs:
    processed_vcf:
        type: File
        outputBinding:
            glob: "add_gt.vcf"

