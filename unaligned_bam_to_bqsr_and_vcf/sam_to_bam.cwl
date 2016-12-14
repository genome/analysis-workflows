#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'convert SAM to BAM'
baseCommand: ["/usr/local/bin/samtools", "view"]
arguments: ["-S", "-b"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/samtools-1.3.1-2:2"
stdout: 'Tagged.bam'
inputs:
    sam:
        type: File
        inputBinding:
            position: 1
outputs:
    bam:
        type: stdout
