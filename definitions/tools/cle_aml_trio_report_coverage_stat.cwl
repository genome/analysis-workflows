#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "cle aml_trio hs_metrics coverage_stat report"
baseCommand: ["/usr/bin/perl", "/usr/local/bin/coverage_stat.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle-aml-trio-reports:v1.0"
    - class: ResourceRequirement
      ramMin: 4000
arguments: [
    $(runtime.outdir)/coverage_stat.out
]
inputs:
    normal_roi_hs_metrics:
        type: File
        inputBinding:
            position: -7 
    normal_summary_hs_metrics:
        type: File[]
        inputBinding:
            position: -6 
    tumor_roi_hs_metrics:
        type: File
        inputBinding:
            position: -5
    tumor_summary_hs_metrics:
        type: File[]
        inputBinding:
            position: -4
    followup_roi_hs_metrics:
        type: File
        inputBinding:
            position: -3 
    followup_summary_hs_metrics:
        type: File[]
        inputBinding:
            position: -2 
outputs:
    coverage_stat:
        type: File
        outputBinding:
            glob: "coverage_stat.out"
