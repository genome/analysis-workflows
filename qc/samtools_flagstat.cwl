#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools flagstat"
baseCommand: ["/usr/local/bin/samtools", "flagstat"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/tagged-alignment:1"
stdout: flagstat.out
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [.bai]
outputs:
    flagstats:
        type: stdout
