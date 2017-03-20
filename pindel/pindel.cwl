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
      ramMin: 16000
inputs:
    tumor_bam:
        type: File
    normal_bam:
        type: File
    tumor_bam_index:
        type: File
    normal_bam_index:
        type: File
    reference:
        type: File
        inputBinding:
            prefix: "-f"
            position: 1
        secondaryFiles: [".fai"]
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
