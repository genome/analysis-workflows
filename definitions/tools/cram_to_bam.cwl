#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools view cram to bam"
baseCommand: ["/usr/local/bin/samtools", "view", "-b"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
arguments:
        ["-o", { valueFrom: $(runtime.outdir)/$(inputs.cram.nameroot).bam }]
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
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
