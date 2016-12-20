#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect alignment summary metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectAlignmentSummaryMetrics"]
arguments:
    ["OUTPUT=", { valueFrom: $(runtime.outdir)/AlignmentSummaryMetrics.txt }]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "INPUT="
    reference:
        type: File
        inputBinding:
            prefix: "REFERENCE_SEQUENCE="
outputs:
    alignment_summary_metrics:
        type: File
        outputBinding:
            glob: "AlignmentSummaryMetrics.txt"
