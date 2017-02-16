#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel v0.2.5b8"
baseCommand: ["/usr/bin/pindel", "-i", "pindel.config"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/pindel-cwl:v0.2.5b8"
    - class: ResourceRequirement
      ramMin: 16000
    - class: InitialWorkDirRequirement
      listing:
          - entryname: tumor.bam
            entry: $(inputs.tumor_bam)
          - entryname: normal.bam
            entry: $(inputs.normal_bam)
          - entryname: tumor.bam.bai
            entry: $(inputs.tumor_bam_index)
          - entryname: normal.bam.bai
            entry: $(inputs.normal_bam_index)
          - entryname: 'pindel.config'
            entry: |
                normal.bam $(inputs.insert_size)    NORMAL
                tumor.bam  $(inputs.insert_size)    TUMOR
arguments:
    ["-w", "20",
     "-o", "all"]
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
