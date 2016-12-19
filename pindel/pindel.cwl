#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel v0.2.5b8"
baseCommand: ["/usr/bin/pindel", "-i", "pindel.config"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/pindel-cwl:v0.2.5b8"
    - class: InitialWorkDirRequirement
      listing:
          - entryname: 'pindel.config'
            entry: |
                $(inputs.normal_bam.path) $(inputs.insert_size)    NORMAL
                $(inputs.tumor_bam.path)  $(inputs.insert_size)    TUMOR
arguments:
    ["-w", "10",
     "-o", "all"]
inputs:
    tumor_bam:
        type: File
        secondaryFiles: .bai
    normal_bam:
        type: File     
        secondaryFiles: .bai
    reference:
        type: File
        inputBinding:
            prefix: "-f"
            position: 1
        secondaryFiles: [".fai"]
    interval_list:
        type: File
        inputBinding:
            prefix: "-j"
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
