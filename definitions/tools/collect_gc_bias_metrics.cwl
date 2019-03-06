#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect gc bias metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectGcBiasMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/GcBiasMetrics.txt },
    "CHART=", { valueFrom: $(runtime.outdir)/GcBiasMetricsChart.pdf },
    "S=", { valueFrom: $(runtime.outdir)/GcBiasMetricsSummary.txt } ]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: mgibio/cle:v1.3.1
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "I="
        secondaryFiles: [^.bai]
    reference:
        type: string
        inputBinding:
            prefix: "R="
    metric_accumulation_level:
        type: string
        inputBinding:
            prefix: "METRIC_ACCUMULATION_LEVEL="
outputs:
    gc_bias_metrics:
        type: File
        outputBinding:
            glob: "GcBiasMetrics.txt"
    gc_bias_metrics_chart:
        type: File
        outputBinding:
            glob: "GcBiasMetricsChart.pdf"
    gc_bias_metrics_summary:
        type: File
        outputBinding:
            glob: "GcBiasMetricsSummary.txt"
