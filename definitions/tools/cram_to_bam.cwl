#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools view cram to bam"
baseCommand: ["/opt/samtools/bin/samtools", "view", "-b"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
arguments:
        ["-o", { valueFrom: $(runtime.outdir)/output.bam }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-T"
            position: 1
    cram:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [^.crai]
outputs:
    bam:
        type: File
        outputBinding:
            glob: "output.bam"
