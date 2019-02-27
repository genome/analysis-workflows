#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run SURVIVOR to merge SV calls"

baseCommand: ["/bin/bash", "/usr/bin/survivor_merge_helper.sh"]

requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/survivor-cwl:1.0.6.1"
    - class: ResourceRequirement
      ramMin: 2000
      coresMin: 1

inputs:
    vcfs:
        type: File[]
        inputBinding:
            position: 1
            itemSeparator: ","
        doc: "Array of VCFs to merge SV calls of"
    max_distance_to_merge:
        type: int
        inputBinding:
            position: 2
        doc: "Maximum distance of variants to consider for merging"
    minimum_sv_calls:
        type: int
        inputBinding:
            position: 3
        doc: "Minimum number of sv calls needed to be merged"
    same_type:
        type: int
        inputBinding:
            position: 4
        doc: "Require merged SVs to be of the same type, 1=yes, 0=no"
    same_strand:
        type: int
        inputBinding:
            position: 5
        doc: "Require merged SVs to be on the same strand, 1=yes, 0=no"
    estimate_sv_distance:
        type: int
        inputBinding:
            position: 6
        doc: "Estimate distance based on the size of SV, 1=yes, 0=no"
    minimum_sv_size:
        type: int
        inputBinding:
            position: 7
        doc: "Minimum size of SVs to merge"
    cohort_name:
        type: string?
        inputBinding:
            position: 8
        default: "SURVIVOR-sv-merged.vcf"
        doc: "Used to generate the output file name"

outputs:
  merged_vcf:
    type: File
    outputBinding:
      glob: "$(inputs.cohort_name)"

