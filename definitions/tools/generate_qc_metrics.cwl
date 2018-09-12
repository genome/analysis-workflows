#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Picard: RNA Seq Metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/opt/picard/picard.jar", "CollectRnaSeqMetrics"]
requirements:
    - class: DockerRequirement
      dockerPull: mgibio/rnaseq
arguments: [ {valueFrom: "O=$(runtime.outdir)/rna_metrics.txt"},
             {valueFrom: "CHART=$(runtime.outdir)/rna_metrics.pdf"} ]
inputs:
    refFlat:
        type: File
        inputBinding:
            prefix: "REF_FLAT="
            separate: false
    ribosomal_intervals:
        type: File
        inputBinding:
            prefix: "RIBOSOMAL_INTERVALS="
            separate: false
    strand:
        type: string
        inputBinding:
            prefix: "STRAND="
            separate: false
    bam:
        type: File
        inputBinding:
            prefix: "I="
            separate: false
outputs:
    metrics:
        type: File
        outputBinding:
            glob: "rna_metrics.txt"
    chart:
        type: File
        outputBinding:
            glob: "rna_metrics.pdf"
