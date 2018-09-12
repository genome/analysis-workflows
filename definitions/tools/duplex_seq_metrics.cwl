#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'collect duplex seq metrics'
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/fgbio-0.5.0.jar", "CollectDuplexSeqMetrics"]
arguments:
    ["--output", { valueFrom: "$(runtime.outdir)/duplex_seq.metrics"} ]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 1000
    - class: DockerRequirement
      dockerPull: mgibio/dna-alignment
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
