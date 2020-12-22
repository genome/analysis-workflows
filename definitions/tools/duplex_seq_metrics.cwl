#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'collect duplex seq metrics'
baseCommand: ["/usr/local/bin/fgbio", "CollectDuplexSeqMetrics"]
arguments:
    ["--output", { valueFrom: "$(runtime.outdir)/duplex_seq.metrics"} ]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 1000
    - class: DockerRequirement
      dockerPull: quay.io/biocontainers/fgbio:1.3.0--0
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "--input"
    intervals:
        type: File?
        inputBinding:
            prefix: "--intervals"
    description:
        type: string
        inputBinding:
            prefix: "--description"
outputs:
    duplex_seq_metrics:
        type: File[]
        outputBinding:
            glob: "duplex_seq.metrics.*"
