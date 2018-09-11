#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'umi alignment'
baseCommand: ["/bin/bash", "/usr/bin/umi_alignment.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 8
    - class: DockerRequirement
      dockerPull: mgibio/dna-alignment
arguments:
    - position: 3
      valueFrom: $(runtime.cores)
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
    reference:
        type: string
        inputBinding:
            position: 2
outputs:
    aligned_bam:
        type: File
        secondaryFiles: [^.bai]
        outputBinding:
            glob: "aligned.bam"
