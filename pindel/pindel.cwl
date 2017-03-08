#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel v0.2.5b8"
arguments: [
    "echo", $(inputs.normal_bam.path), "\t", $(inputs.insert_size), "\tNORMAL", { valueFrom: " > ", shellQuote: false }, "pindel.config",
    { valueFrom: " && ", shellQuote: false },
    "echo", $(inputs.tumor_bam.path), "\t", $(inputs.insert_size), "\tTUMOR", { valueFrom: " >> ", shellQuote: false }, "pindel.config",
    { valueFrom: " && ", shellQuote: false },
    "/usr/bin/pindel", "-i", "pindel.config", "-w", "20", "-o", "all"
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
