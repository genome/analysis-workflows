#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools view cram to bam"
baseCommand: ["/opt/samtools/bin/samtools", "view", "-b"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
arguments:
        ["-o", { valueFrom: $(runtime.outdir)/$(inputs.cram.nameroot).bam }]
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
            glob: "$(inputs.cram.nameroot).bam"
