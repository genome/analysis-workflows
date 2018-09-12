#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel somatic filter v1"
arguments: [
    "/usr/bin/perl", "/usr/bin/write_pindel_filter_config.pl", $(inputs.pindel_output_summary.path), $(inputs.reference), $(runtime.outdir),
    { valueFrom: " && ", shellQuote: false },
    "/usr/bin/perl", "/usr/bin/somatic_indelfilter.pl", "filter.config"
]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: mgibio/cle
inputs:
    reference:
        type: string
    pindel_output_summary:
        type: File
outputs:
    vcf:
        type: File
        outputBinding:
            glob: "pindel.out.vcf"
