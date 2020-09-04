#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools flagstat"
baseCommand: ["/opt/samtools/bin/samtools", "flagstat"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
stdout: "$(inputs.cram.basename).flagstat"
inputs:
    reference:
        type: File
        secondaryFiles: [.fai]
        inputBinding:
            position: 1
            prefix: '--reference'
    cram:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [^.crai]
outputs:
    flagstats:
        type: stdout
