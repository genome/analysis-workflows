#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Picard: RNA Seq Metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/opt/picard/picard.jar", "CollectRnaSeqMetrics"]
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
    metrics_output_file:
        type: string
        inputBinding:
            prefix: "O="
            separate: false
    chart_output_file:
        type: string
        inputBinding:
            prefix: "CHART="
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
            glob: $(inputs.metrics_output_file)
    chart:
        type: File
        outputBinding:
            glob: $(inputs.chart_output_file)
