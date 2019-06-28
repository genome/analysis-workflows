#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "cle aml_trio full variants report"
baseCommand: ["/usr/bin/perl", "/usr/local/bin/full_variant_report.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle-aml-trio-reports:v1.0"
    - class: ResourceRequirement
      ramMin: 4000
arguments: [
    $(runtime.outdir)/full_variant_report.out
]
inputs:
    variant_tsv:
        type: File
        inputBinding:
            position: -5
    followup_snv_bam_readcount:
        type: File
        inputBinding:
            position: -4
    followup_indel_bam_readcount:
        type: File
        inputBinding:
            position: -3
    pindel_region_vcf:
        type: File
        inputBinding:
            position: -2
outputs:
    full_variant_report:
        type: File
        outputBinding:
            glob: "full_variant_report.out"

