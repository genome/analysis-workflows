#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "varscan v2.4.2 somatic"
baseCommand: "/usr/bin/varscan_helper.sh"
requirements:
    - class: ResourceRequirement
      ramMin: 12000
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
inputs:
    tumor_cram:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [^.crai]
    normal_cram:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [^.crai]
    reference:
        type: string
        inputBinding:
            position: 3
    strand_filter:
        type: int?
        default: 0
        inputBinding:
            position: 4
    min_coverage:
        type: int?
        default: 8
        inputBinding:
            position: 5
    min_var_freq:
        type: float?
        default: 0.1
        inputBinding:
            position: 6
    p_value:
        type: float?
        default: 0.99
        inputBinding:
            position: 7
    roi_bed:
        type: File?
        inputBinding:
            position: 8
outputs:
    snvs:
        type: File
        outputBinding:
            glob: "output.snp.vcf"
    indels:
        type: File
        outputBinding:
            glob: "output.indel.vcf"
