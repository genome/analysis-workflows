#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools flagstat"
baseCommand: ["/opt/samtools/bin/samtools", "flagstat"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
stdout: flagstat.out
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [^.bai]
outputs:
    flagstats:
        type: stdout
