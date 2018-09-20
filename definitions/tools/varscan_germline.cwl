#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "varscan v2.4.2 germline"
baseCommand: "/usr/bin/varscan_germline_helper.sh"
requirements:
    - class: ResourceRequirement
      ramMin: 12000
    - class: DockerRequirement
      dockerPull: "mgibio/varscan_helper-cwl:1.0.0"
arguments:
    - position: 9
      valueFrom: "$(runtime.outdir)/output.vcf"
inputs:
    cram:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [^.crai]
    reference:
        type: string
        inputBinding:
            position: 2
    strand_filter:
        type: int?
        default: 0
        inputBinding:
            position: 3
    min_coverage:
        type: int?
        default: 8
        inputBinding:
            position: 4
    min_var_freq:
        type: float?
        default: 0.1
        inputBinding:
            position: 5
    min_reads:
        type: int?
        default: 2
        inputBinding:
            position: 6
    p_value:
        type: float?
        default: 0.99
        inputBinding:
            position: 7
    sample_name:
        type: string
        inputBinding:
            position: 8
    roi_bed:
        type: File?
        inputBinding:
            position: 10
outputs:
    variants:
        type: File
        outputBinding:
            glob: "output.vcf"
