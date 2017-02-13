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
      coresMin: 4
    - class: InitialWorkDirRequirement
      listing:
          - entryname: 'pindel.config'
            entry: |
                $(inputs.normal_bam.path) $(inputs.insert_size)    NORMAL
                $(inputs.tumor_bam.path)  $(inputs.insert_size)    TUMOR
arguments:
    ["-w", "20",
     "-T", "4",
     "-o", "all"]
inputs:
    tumor_bam:
        type: File
        secondaryFiles: ["^.bai"]
    normal_bam:
        type: File
        secondaryFiles: ["^.bai"]
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
