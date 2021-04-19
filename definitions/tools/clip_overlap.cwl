#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'clip overlapping reads'
baseCommand: ["/usr/local/bin/fgbio", "ClipBam"]
arguments:
    ["--clipping-mode", "Hard", "--clip-overlapping-reads", "true",
    "--output", { valueFrom: "$(runtime.outdir)/clipped.bam"} ]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: quay.io/biocontainers/fgbio:1.3.0--0
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "--input"
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "--ref"
outputs:
    clipped_bam:
        type: File
        secondaryFiles: [^.bai]
        outputBinding:
            glob: "clipped.bam"
