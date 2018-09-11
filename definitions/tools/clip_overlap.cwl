#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'clip overlapping reads'
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/fgbio-0.5.0.jar", "ClipBam"]
arguments:
    ["--soft-clip", "false", "--clip-overlapping-reads", "true",
    "--output", { valueFrom: "$(runtime.outdir)/clipped.bam"} ]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/dna-alignment
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "--input"
    reference:
        type: string
        inputBinding:
            prefix: "--ref"
outputs:
    clipped_bam:
        type: File
        outputBinding:
            glob: "clipped.bam"
