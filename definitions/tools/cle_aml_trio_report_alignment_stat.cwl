#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "cle aml_trio alignment_stat report"
baseCommand: ["/usr/bin/perl", "/usr/local/bin/alignment_stat.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle-aml-trio-reports:v1.0"
    - class: ResourceRequirement
      ramMin: 4000
arguments: [
    $(runtime.outdir)/alignment_stat.out
]
inputs:
    normal_alignment_summary_metrics:
        type: File
        inputBinding:
            position: -4
    tumor_alignment_summary_metrics:
        type: File
        inputBinding:
            position: -3
    followup_alignment_summary_metrics:
        type: File
        inputBinding:
            position: -2
outputs:
    alignment_stat:
        type: File
        outputBinding:
            glob: "alignment_stat.out"
