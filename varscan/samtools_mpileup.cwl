#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools 1.3.1 somatic mpileup"
baseCommand: "mpileup"
arguments:
    ["--no-baq"]
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "-f"
            position: 2
        secondaryFiles: [.fai]
    normal_bam:
        type: File
        inputBinding:
            position: 3
    tumor_bam:
        type: File
        inputBinding:
            position: 4
outputs:
    mpileup:
        type: stdout
