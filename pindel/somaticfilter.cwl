#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel somatic filter v1"
arguments: [
    "/usr/bin/perl", "write_pindel_filter_config.pl", $(inputs.pindel_output_summary.path), $(inputs.reference.path), $(runtime.outdir)",
    { valueFrom: " && ", shellQuote: false },
    "/usr/bin/perl", "/usr/bin/somatic_indelfilter.pl", "filter.config"
]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
    - class: ShellCommandRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [".fai"]
    pindel_output_summary:
        type: File
outputs:
    vcf:
        type: File
        outputBinding:
            glob: "pindel.out.vcf"
