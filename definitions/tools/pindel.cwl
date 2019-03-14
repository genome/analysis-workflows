#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel v0.2.5b8"
arguments: [
    "/usr/bin/perl", "/usr/bin/pindel_helper.pl",
    $(inputs.normal_bam.path), $(inputs.tumor_bam.path), $(inputs.insert_size)
]
requirements:
    - class: ResourceRequirement
      ramMin: 64000
      tmpdirMin: 100000
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.3.1"
inputs:
    tumor_bam:
        type: File
        secondaryFiles: ["^.bai"]
    normal_bam:
        type: File
        secondaryFiles: ["^.bai"]
    reference:
        type: string
        inputBinding:
            prefix: "-f"
            position: 1
    chromosome:
        type: string
        inputBinding:
            prefix: "-c"
            position: 2
    insert_size:
        type: int
        default: 400
outputs:
    deletions:
        type: File
        outputBinding:
            glob: "all_D"
    insertions:
        type: File
        outputBinding:
            glob: "all_SI"
    tandems:
        type: File
        outputBinding:
            glob: "all_TD"
    long_insertions:
        type: File
        outputBinding:
            glob: "all_LI"
    inversions:
        type: File
        outputBinding:
            glob: "all_INV"
