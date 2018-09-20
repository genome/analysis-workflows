#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect alignment summary metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectAlignmentSummaryMetrics"]
arguments:
    ["OUTPUT=", { valueFrom: $(runtime.outdir)/AlignmentSummaryMetrics.txt }]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
    - class: DockerRequirement
      dockerPull: "mgibio/picard-cwl:2.18.1"
inputs:
    cram:
        type: File
        inputBinding:
            prefix: "INPUT="
        secondaryFiles: [^.crai]
    reference:
        type: string
        inputBinding:
            prefix: "REFERENCE_SEQUENCE="
    metric_accumulation_level:
        type: string
        inputBinding:
            prefix: "METRIC_ACCUMULATION_LEVEL="
outputs:
    alignment_summary_metrics:
        type: File
        outputBinding:
            glob: "AlignmentSummaryMetrics.txt"
