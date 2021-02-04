#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools flagstat"
baseCommand: ["/usr/local/bin/samtools", "flagstat"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
stdout: "$(inputs.bam.basename).flagstat"
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [^.bai]
outputs:
    flagstats:
        type: stdout
