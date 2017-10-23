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
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
    - class: ResourceRequirement
      ramMin: 16000
inputs:
    cram:
        type: File
        inputBinding:
            prefix: "I="
        secondaryFiles: [^.crai]
    reference:
        type: string
        inputBinding:
            prefix: "R="
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
