#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'umi realignment'
baseCommand: ["/bin/bash", "/usr/bin/umi_realignment.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 8
    - class: DockerRequirement
      dockerPull: "mgibio/dna-alignment:1.0.0"
arguments:
    - position: 3
      valueFrom: $(runtime.cores)
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
    reference:
        type:
            - string
            - File
        secondaryFiles: [.amb, .ann, .bwt, .pac, .sa]
        inputBinding:
            position: 2
outputs:
    consensus_aligned_bam:
        type: File
        secondaryFiles: [^.bai]
        outputBinding:
            glob: "realigned.bam"
